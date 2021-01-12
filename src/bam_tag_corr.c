#include "utils.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "bam_pool.h"
#include "htslib/kstring.h"
#include "number.h"
#include "htslib/thread_pool.h"
#include "thread.h"
#include "dict.h"
#include "htslib/khash.h"
#include "dna_pool.h"
#include "multi_index_list.h"

// maximal hamming distance of two similar UMI
#define UMI_E 1

static struct args {
    const char  * input_fname;
    const char  * output_fname;

    const char  * tag;
    int           n_block;
    char        **blocks;
    struct dict **blkidx;

    struct index_list *idxlst;
    
    const char  * new_tag;

    int           cr_method; // correct umi like cellrange
    
    htsFile     * in;
    htsFile     * out;

    bam_hdr_t   * hdr;
    
    int           n_thread;
    int           file_th;
    
    int           chunk_size;
    
    int           dist;
} args = {
    .input_fname  = NULL,
    .output_fname = NULL,
    .tag          = NULL,
    .new_tag      = NULL,
    .n_block      = 0,
    .blocks       = NULL,
    .blkidx       = NULL,
    .idxlst       = NULL,
    .cr_method    = 0,
    .in           = NULL,
    .out          = NULL,
    .hdr          = NULL,
    .n_thread     = 5,
    .file_th      = 4,

    .chunk_size   = 1000000, //1M

    .dist         = UMI_E,
};

static void memory_release()
{
    bam_hdr_destroy(args.hdr);
    sam_close(args.in);
    sam_close(args.out);
    int i;
    for (i = 0; i < args.idxlst->l; ++i) {
        struct PISA_dna_pool *p = args.idxlst->mi[i].data;
        PISA_dna_destroy(p);
    }
    index_list_destroy(args.idxlst);

    for (i = 0; i < args.n_block; ++i) {
        free(args.blocks[i]);
        dict_destroy(args.blkidx[i]);
    }
    free(args.blocks);
    free(args.blkidx);
}

int *sam_tag_values(bam1_t *b, int n,  char **blocks)
{
    int *idx = malloc(n*sizeof(int));
    int i;
    for (i = 0; i < n; ++i) {
        const char *v = (const char*)bam_aux_get(b, blocks[i]);
        if (!v) {
            free(idx); return NULL;
        }
        int ret;
        ret = dict_query(args.blkidx[i], v+1);
        if (ret < 0) ret = dict_push(args.blkidx[i], v+1);
        idx[i] = ret;
        assert(idx[i]>=0);
    }

    return idx;
}


void build_index()
{
    LOG_print("Building index ..");
    double t_real;
    t_real = realtime();
    
    //struct dict *cell_bc = dict_init();
    //dict_set_value(cell_bc);
    
    htsFile *fp = hts_open(args.input_fname, "r");
    
    CHECK_EMPTY(fp, "%s : %s.", args.input_fname, strerror(errno));
    htsFormat type = *hts_get_format(fp);
    if (type.format != bam && type.format != sam) error("Unsupported input format, only support BAM/SAM/CRAM format.");
    
    bam_hdr_t *hdr = sam_hdr_read(fp);
    hts_set_threads(fp, args.file_th);
    bam1_t *b = bam_init1();
    int i;
    for (;;) {
        if (sam_read1(fp, hdr, b) < 0) break;
        bam1_core_t *c = &b->core;
        if (c->flag & BAM_FQCFAIL ||
            c->flag & BAM_FSECONDARY ||
            c->flag & BAM_FSUPPLEMENTARY ||
            c->flag & BAM_FUNMAP ||
            c->flag & BAM_FDUP) continue;
        
        const char *umi = (const char *)bam_aux_get(b, args.tag);
        if (!umi) continue;
        int *idx = sam_tag_values(b, args.n_block, args.blocks);
        if (idx == NULL) continue;
        
        struct multi_index *mi = index_query(args.idxlst, args.n_block, idx);
        struct PISA_dna_pool *u;
        if (mi == NULL) {
            u = PISA_dna_pool_init();
            index_push1(args.idxlst, args.n_block, idx, u);
        }
        else {
            u = (struct PISA_dna_pool*) mi->data;
            free(idx);
        }
        assert(u);
        PISA_dna_push(u, umi+1);
    }
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    /*
    for (i = 0; i < args.idxlst->l; ++i) {
        struct multi_index *mi = &args.idxlst->mi[i];
        struct PISA_dna_pool *u = mi->data;
        int j, k;
        for (k = 0; k <u->l; ++k) {
            for (j = 0; j < args.n_block; ++j) {
                fprintf(stdout, "%s\t", dict_name(args.blkidx[j], mi->mi[j]));
            }
            char *str = dna_decode_str(&u->data[k],u->len/2);
            fprintf(stdout, "%s\t%d\n", str, u->data[k].count);
            free(str);
        }
        
    }
    */
    // correction
    // todo: multi-threads
    for (i = 0; i < args.idxlst->l; ++i) {
        struct PISA_dna_pool *u = args.idxlst->mi[i].data;
        assert(u);
        dna_pool_corr(u, args.dist);
    }
    /*
    for (i = 0; i < args.idxlst->l; ++i) {
        struct multi_index *mi = &args.idxlst->mi[i];
        struct PISA_dna_pool *u = mi->data;
        int j, k;
        for (k = 0; k <u->l; ++k) {
            fprintf(stdout, "corr\t");
            for (j = 0; j < args.n_block; ++j) {
                fprintf(stdout, "%s\t", dict_name(args.blkidx[j], mi->mi[j]));
            }
            char *str = dna_decode_str(&u->data[k],u->len/2);
            fprintf(stdout, "%s\t%d\n", str, u->data[k].count);
            free(str);
        }
        
    }
    */
    if (args.cr_method) {
    }
        
    /*
    //struct tpool *tp = tpool_init(args.n_thread, args.n_thread*2, 1);

    for (i = 0; i < dict_size(cell_bc); ++i) {
        struct bc_corr *bc0 = dict_query_value(cell_bc, i);
        build_index1(bc0);
        //tpool_add_work(tp, build_index1, bc0);
    }
    //tpool_destroy(tp);
    */
    
    LOG_print("Build time : %.3f sec", realtime() - t_real);
}

static int parse_args(int argc, char **argv)
{
    if (argc ==1) return 1;
    
    const char *block_tags = NULL;
    const char *file_th = NULL;
    const char *thread = NULL;
    const char *distance = NULL;
    
    int i;
    for (i = 1; i < argc;) {
        const char *a = argv[i++];
        const char **var = 0;

        if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-h") == 0) return 1;
        else if (strcmp(a, "-tag") == 0) var = &args.tag;
        else if (strcmp(a, "-tags-block") == 0) var = &block_tags;
        else if (strcmp(a, "-@") == 0) var = &file_th;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-new-tag") == 0) var = &args.new_tag;
        else if (strcmp(a, "-e") == 0) var = &distance;
        else if (strcmp(a, "-cr") == 0) {
            args.cr_method = 1;
            continue;
        }

        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (a[0] == '-' && a[1] != '\0') error("Unknown option, %s", a);
        
        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument : %s",a);
    }

    CHECK_EMPTY(args.output_fname, "-o need to be set.");
    CHECK_EMPTY(args.input_fname, "No input bam.");
    CHECK_EMPTY(args.tag, "-tag need to be set.");
    CHECK_EMPTY(block_tags, "-tags-block need to be set.");
    
    kstring_t str = {0,0,0};
    kputs(block_tags, &str);
    int *s = ksplit(&str, ',', &args.n_block);
    assert(args.n_block >0);
    args.blocks = malloc(args.n_block*sizeof(char*));
    args.blkidx = malloc(args.n_block*sizeof(void*));
    for (i = 0; i < args.n_block; ++i) {
        args.blocks[i] = strdup(str.s+s[i]);
        args.blkidx[i] = dict_init();
    }
    free(s);
    free(str.s);

    if (file_th) args.file_th = str2int((char*)file_th);
    if (thread) args.n_thread = str2int((char*)thread);
    if (distance) args.dist = str2int((char*)distance);
    if (args.dist < 1) error("Hamming distance of similar barcodes greater than 0 is required.");

    args.idxlst = index_list_init();
    
    build_index();
    return 0;
}

int update_new_tag(bam1_t *b)
{
    char *umi = (char*)bam_aux_get(b, args.tag);
    if (!umi) return 0;

    int *idx = sam_tag_values(b, args.n_block, args.blocks);
    if (idx == NULL) return 0;

    struct multi_index *mi = index_query(args.idxlst, args.n_block, idx);
    assert(mi && mi->data);

    struct PISA_dna_pool *u = mi->data;

    char *new_umi = dna_corr(u, umi+1);    
    if (!new_umi) {
        int i;
        for (i = 0; i < args.n_block; ++i) {
            debug_print("idx : %d", idx[i]);
        }
        error();
        return 0;
    }
    free(idx);    
    if (args.new_tag)
        bam_aux_append(b, args.new_tag, 'Z', strlen(new_umi)+1, (uint8_t*)new_umi);
    else
        memcpy(umi+1, new_umi, strlen(new_umi)); // since it is equal length, just reset the memory..
    
    free(new_umi);
    return 1;
}

static void *run_it(void *data)
{
    struct bam_pool *p = (struct bam_pool*)data;
    int i;
    int c = 0;
    for (i = 0; i < p->n; ++i) {
        bam1_t *b = &p->bam[i];
        c += update_new_tag(b);
    }
    
    return p;
}
static void write_out(struct bam_pool *p)
{
    int i;
    for (i = 0; i < p->n; ++i)        
        if (sam_write1(args.out, args.hdr, &p->bam[i]) == -1) error("Failed to write SAM.");
    bam_pool_destory(p);
}

extern int bam_corr_usage();

int bam_corr_umi(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return bam_corr_usage();

    args.in  = hts_open(args.input_fname, "r");
    args.hdr = sam_hdr_read(args.in);
    CHECK_EMPTY(args.hdr, "Failed to open header.");
    
    args.out = hts_open(args.output_fname, "bw");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));
    
    if (sam_hdr_write(args.out, args.hdr)) error("Failed to write SAM header.");
    
    hts_set_threads(args.out, args.file_th); // write file in multi-threads

    /*
    for (;;) {
        struct bam_pool *b = bam_pool_create();
        bam_read_pool(b, args.in, args.hdr, args.chunk_size);
            
        if (b == NULL) break;
        if (b->n == 0) { free(b->bam); free(b); break; }
        
        b = run_it(b);
        write_out(b);   
    }
    */

    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;

    for (;;) {
        struct bam_pool *b = bam_pool_create();
        bam_read_pool(b, args.in, args.hdr, args.chunk_size);
            
        if (b == NULL) break;
        if (b->n == 0) { free(b->bam); free(b); break; }
        
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b, 1);
            if ((r = hts_tpool_next_result(q))) {
                struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
                write_out(d);   
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);        
    }
    
    hts_tpool_process_flush(q);

    while ((r = hts_tpool_next_result(q))) {
        struct bam_pool *d = (struct bam_pool*)hts_tpool_result_data(r);
        write_out(d);
        hts_tpool_delete_result(r, 0);
    }
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    memory_release();
    LOG_print("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - t_real, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    //LOG_print("%d records updated.", args.update_count);
    return 0;
}
    
#ifdef CORR_UMI
int main(int argc, char **argv)
{
    return bam_corr_umi(argc, argv);
}

#endif
