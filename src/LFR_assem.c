#include "utils.h"
#include "rld0.h"
#include "rle.h"
#include "mrope.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"
#include "fastq.h"
#include "ksw.h"
#include "fml.h"
#include "number.h"
#include "mag.h"

KSEQ_INIT(gzFile, gzread)

extern int ksa_sa(const unsigned char *T, int *SA, int n, int k);
extern int ksa_bwt(unsigned char *T, int n, int k);

static unsigned char seq_nt6_table[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

static void revcomp6(uint8_t *s, int l)
{
	int i;
	for (i = 0; i < l>>1; ++i) {
		int tmp = s[l-1-i];
		tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp;
		s[l-1-i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
		s[i] = tmp;
	}
	if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i];
}
static uint8_t *enc_str(char *s, int l)
{
    uint8_t *e = malloc(l+1);
    int i;
    for (i = 0; i < l; ++i) e[i] = seq_nt6_table[(int)s[i]];
    e[l] = 0;
    return e;
}
static void assem_opt_init(fml_opt_t *opt)
{
    fml_opt_init(opt);
    opt->min_asm_ovlp = 10;
    opt->ec_k = -1;
    opt->mag_opt.flag = MAG_F_NO_SIMPL | MAG_F_AGGRESSIVE;
};
static struct args {
    const char *input_fname;
    const char *output_fname;
    int n_tag;
    char **tags;    
    int pair;
    int min_ovlp;
    int n_thread;
    char *last_name;

    gzFile fp_in;
    FILE *out;
    kseq_t *ks;
    fml_opt_t *assem_opt;
    int l_seed;
    uint8_t *seed;
    uint64_t filter_block;
    uint64_t assem_block;    
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .n_tag = 0,
    .tags = NULL,
    .pair = 0,
    .min_ovlp = 10,
    .n_thread = 1,
    .last_name = NULL,
   
    .fp_in = NULL,
    .out = NULL,
    .ks = NULL,
    .assem_opt = NULL,
    .l_seed = 0,
    .seed = NULL,
    .filter_block = 0,
    .assem_block = 0,
};

static int usage()
{
    fprintf(stderr, "scLFR_assem in.fq\n");
    fprintf(stderr, "   -t         Threads.\n");
    fprintf(stderr, "   -o         Output fastq.\n");
    fprintf(stderr, "   -tag       Tags of read block.\n");
    fprintf(stderr, "   -ss        Seed sequence for scLFR library. Usually be fixed adaptor sequences.\n");
    // fprintf(stderr, "   -p         Input fastq is smart paired.\n");
    return 1;
}
static void memory_release()
{
    gzclose(args.fp_in);
    kseq_destroy(args.ks);
    fclose(args.out);    
}
static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;

    int i;
    const char *thread = NULL;
    const char *tags = NULL;
    const char *seed = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        const char **var = 0;
        if (strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0) return 1;
        else if (strcmp(a, "-t") == 0) var = &thread;
        else if (strcmp(a, "-tag") == 0) var = &tags;
        else if (strcmp(a, "-o") == 0) var = &args.output_fname;
        else if (strcmp(a, "-ss") == 0) var = &seed;
        else if (strcmp(a, "-p") == 0) {
            args.pair = 1;
            continue;
        }

        if (var != 0) {
            *var = argv[i++];
            continue;
        }

        if (args.input_fname == NULL) {
            args.input_fname = a;
            continue;
        }

        error("Unknown argument : %s", a);
    }
    if (args.input_fname == NULL ) error("No input fastq specified.");
    if (tags == NULL) error("-tag must be set.");
    if (seed) { //error("Seed sequence for scLFR must be set with -ss.");
        args.l_seed = strlen(seed);
        if (args.l_seed < 6) error("Seed is too short, need >= 6bp.");
        args.seed = enc_str(seed, args.l_seed);
    }
    
    kstring_t str = {0,0,0};
    kputs(tags, &str);
    int n;
    int *s = ksplit(&str, ',', &n);
    args.n_tag = n;
    args.tags = malloc(n*sizeof(char*));
    for (i = 0; i <n; ++i) args.tags[i] = strdup(str.s+s[i]);
    free(s); free(str.s);

    if (thread) args.n_thread = str2int((char*)thread);
    if (args.n_thread < 1) args.n_thread = 1;

    args.out = args.output_fname == NULL ? stdout : fopen(args.output_fname, "w");
    CHECK_EMPTY(args.out, "%s : %s.", args.output_fname, strerror(errno));

    args.fp_in = gzopen(args.input_fname, "r");
    CHECK_EMPTY(args.fp_in, "%s : %s.", args.input_fname, strerror(errno));

    args.ks = kseq_init(args.fp_in);
    args.assem_opt = malloc(sizeof(fml_opt_t));
    assem_opt_init(args.assem_opt);
    return 0;
}

struct read_block {
    char *name;
    bseq1_t *b;
    int n, m;
};

static void read_block_destory(struct read_block *r)
{
    if (r->name) free(r->name);
    int i;
    for (i = 0; i < r->n; ++i) {
        free(r->b[i].seq);
        if (r->b[i].qual) free(r->b[i].qual);
    }
    free(r->b);
    free(r);    
}
static char *generate_names(char **names)
{
    kstring_t str = {0,0,0};
    int i;
    for (i = 0; i < args.n_tag; ++i) kputs(names[i], &str);
    for (i = 0; i < args.n_tag; ++i) {
        kputs("|||", &str);
        kputs(args.tags[i], &str);
        kputs("|||", &str);
        kputs(names[i], &str);
    }
    return str.s;
}
static struct read_block *read_block()
{
    struct read_block *b =  malloc(sizeof(*b));
    memset(b, 0, sizeof(*b));
    b->m = 4;
    b->b = malloc(sizeof(bseq1_t)*b->m);

    if (args.last_name) {
        b->name = strdup(args.last_name);

        bseq1_t *bb = &b->b[0];
        memset(bb, 0, sizeof(*bb));
        bb->seq = strdup(args.ks->seq.s);
        if (args.ks->qual.s) bb->qual =  strdup(args.ks->qual.s);
        bb->l_seq = args.ks->seq.l;
        b->n++;
        free(args.last_name);
        args.last_name = NULL;
    }

    // read name:  TAG1TAGS|||TAG1|||VAL1|||TAG2|||VAL2
    while (kseq_read(args.ks)>=0) {
        if (b->n == b->m) {
            b->m = b->m<<1;
            b->b = realloc(b->b, sizeof(bseq1_t)*b->m);
        }

        char **name = fastq_name_pick_tags(args.ks->name.s, args.n_tag, args.tags);
        char *n = generate_names(name);
        int i;
        for (i = 0; i <args.n_tag; ++i) free(name[i]);
        free(name);
        if (b->name == NULL) b->name = strdup(n);
        // if (args.last_name == NULL) args.last_name = strdup(n); // first line
        if (strcmp(n, b->name) == 0) {
            bseq1_t *b1 = &b->b[b->n++];
            memset(b1, 0, sizeof(*b1));
            b1->seq = strdup(args.ks->seq.s);
            if (args.ks->qual.l) b1->qual = strdup(args.ks->qual.s);
            b1->l_seq = args.ks->seq.l;
            free(n);
        }
        else {
            //if (args.last_name) free(args.last_name);
            args.last_name = n;
            break;
        }
    }
    if (b->n == 0) {
        free(b); return NULL;
    }
    return b;
}

struct base_v {
    int l;
    uint8_t *v;
};
void push_str_base(struct base_v *v, char *s)
{
    if (s == NULL) return;
    int l;
    l = strlen(s);
    v->v = realloc(v->v, v->l+(l+1)*2);
    uint8_t *e = enc_str(s, l);
    memcpy(v->v+v->l, e, l+1);
    revcomp6(e, l);
    memcpy(v->v + v->l + l+1, e, l+1);
    v->l = v->l + (l+1)*2;
    free(e);
}

extern int merge_paired(const int l_seq1, const int l_seq2, char const *s1, char const *s2, char const *q1, char const *q2, char **_seq, char **_qual);


static struct base_v *rend_bseq(struct read_block *b)
{
    struct base_v *v = malloc(sizeof(*v));
    v->l = 0; v->v = NULL;

    int i;
    for (i = 0; i < b->n; ) {
        bseq1_t *b1 = &b->b[i++];
        /*
        if (pair) {
            bseq1_t *b2 = &b->b[i++];
            assert(i < b->n);
            char *seq = NULL;
            char *qual = NULL;
            if (PEmerge(b1->l_seq, b2->l_seq, b1->seq, b2->seq, b1->qual, b2->qual, &seq, &qual) == 0) {
                push_str_base(v, seq);
                free(seq);
                if (qual) free(qual);
            }
            else {
                push_str_base(v, b1->seq);
                push_str_base(v, b2->seq);
            }
        }
        else 
        */
        push_str_base(v, b1->seq);
    }    
    return v;
}
/*
static char *rend_utg(char *name, fml_utg_t *utg, int n)
{
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < n; ++i) {
        fml_utg_t *u = &utg[i];
        kputc('@', &str);
        kputw(i, &str);kputc('_', &str);
        kputs(name, &str);
        kputc('\t', &str);
        kputs("assembled",&str);
        kputc('\n', &str);
        kputsn(u->seq, u->len, &str); kputc('\n', &str);
        kputc('+', &str);  kputc('\n', &str);
        kputsn(u->cov, u->len, &str); kputc('\n', &str);
    }
    kputs("", &str);
    return str.s;
}
*/
/*
mrope_t *ropebwt_build(struct fastq_handler *fastq)
{
    
    kstring_t str = {0,0,0};
    
    
}
*/
static int bwt_gen(int asize, int64_t l, uint8_t *s)
{
    if (l <= INT32_MAX) return ksa_bwt(s, l, asize);
    else error("Do not support so many reads.");
}

static rld_t *bwt_enc(int asize, int sbits, int64_t l, const uint8_t *s)
{
	int c;
	int64_t i, k;
	rlditr_t itr;
	rld_t *e;

	e = rld_init(asize, sbits);
	rld_itr_init(e, &itr, 0);
	k = 1; c = s[0];
	for (i = 1; i < l; ++i) {
		if (s[i] != c) {
			rld_enc(e, &itr, k, c);
			c = s[i];
			k = 1;
		} else ++k;
	}
	rld_enc(e, &itr, k, c);
	rld_enc_finish(e, &itr);
	return e;
}
static int check_polyTs(const rld_t *e, int len)
{
    uint64_t k,l, ok, ol;
    int i;
    k = e->cnt[4];
    l = e->cnt[5]-1;
    for (i = len -2; i>=0; --i) {
        rld_rank21(e, k, l, 4, &ok, &ol);
        k = e->cnt[4]+ok;
        l = e->cnt[4]+ol;
        if (k >= l) break;
    }
    if ( k >= l ) return 1;
    return 0;
}
static uint64_t bwt_backward_search(const rld_t *e, int len, const uint8_t *str, uint64_t *sa_beg, uint64_t *sa_end)
{
	uint64_t k, l, ok, ol;
	int i, c;
	c = str[len - 1];
	k = e->cnt[c];
        l = e->cnt[c + 1]-1;
	for (i = len - 2; i >= 0; --i) {
            c = str[i];
            //putchar("$ACGTN"[c]);
            rld_rank21(e, k, l, c, &ok, &ol);
            k = e->cnt[c] + ok;
            l = e->cnt[c] + ol;
            // debug_print("%d\t%c\t%d\t%d",i,"$ACGTN"[c], k, l);
            if (k == l) break;
	}
	if (k == l) return 0;
	*sa_beg = k; *sa_end = l;
	return l - k;
}
static int bwt_ext(const rld_t *e, uint64_t st, uint64_t ed, uint64_t *k0)
{
    uint64_t i;
    int j = 0;
    for (i = st; i <=ed; ++i) {
        uint64_t k = i, ok[6];
        while (1) {
            int c = rld_rank1a(e, k, ok);
            // putchar("$ACGTN"[c]);
            k = e->cnt[c] + ok[c];
            // debug_print("%d\t%c\t%d",i,"$ACGTN"[c], k);
            if (c==0) break;
        }
        k0[j++] = k-1;
        // putchar('\n');
    }
    return j;
}
static rld_t *bwt_build_core(rld_t *e0, int asize, int sbits, int64_t l, uint8_t *s)
{
	rld_t *e;
        bwt_gen(asize, l, s);
        // putchar('\n');                
        e = bwt_enc(asize, sbits, l, s);
	return e;
}

rld_t *bwt_build(struct base_v *v)
{
    return bwt_build_core(0, 6, 3, v->l, v->v);
}
int64_t bwt_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
    uint64_t k = x, *ok;
    ok = alloca(8 * e->asize);
    s->l = 0;
    while (1) {
        int c = rld_rank1a(e, k, ok);
        k = e->cnt[c] + ok[c] -1;
        if (c == 0) return k;
        kputc(c, s);
    }
}
extern int unitig1(aux_t *a, int64_t seed, kstring_t *s, kstring_t *cov, uint64_t end[2], ku128_v nei[2], int *n_reads);

static char *rend_utg(char *name, fml_utg_t *utg, int n)
{
    int i;
    kstring_t str = {0,0,0};
    for (i = 0; i < n; ++i) {
        fml_utg_t *u = &utg[i];
        kputc('>', &str);
        kputw(i, &str);kputc('_', &str);
        kputs(name, &str);
        kputc('\n', &str);
        kputsn(u->seq, u->len, &str); kputc('\n', &str);
    }
    kputs("", &str);
    return str.s;
}

static void *run_it(void *_d)
{
    struct read_block *b = (struct read_block*)_d;
    if (b->n == 0) return NULL;

    // Step 1: Pure clean up of ME sequences, merge short reads 
    struct base_v *v = rend_bseq(b);    

    // Step 2: build eBWT
    rld_t *e = bwt_build(v);
    free(v->v); free(v);
    //debug_print("%s", b->name);
    // Step 3: adaptors and polyTs, if no just skip this block
    if (check_polyTs(e, 10)) goto empty_block;
    if (args.l_seed) {
        uint64_t n;
        uint64_t st = 0, ed = 0;
        n = bwt_backward_search(e, args.l_seed, args.seed, &st, &ed);
        if (n == 0) goto empty_block;
    }
    // Step 4: construct unitigs
    /*
    aux_t a;
    memset(&a, 0, sizeof(a));
    a.e = e;
    a.used = (uint64_t*)calloc((e->mcnt[1] + 63)/64, 8);
    a.bend = (uint64_t*)calloc((e->mcnt[1] + 63)/64, 8);
    a.min_match = args.min_ovlp;
    a.min_merge_len = 0;
    int i, j = 0;
    kstring_t str, cov;
    magv_t z;
    for (i = 1; i < e->mcnt[1]; i+=2) {
        str.l = 0;
        cov.l = 0;
        if (unitig1(&a, i, &str, &cov, z.k, z.nei, &z.nsr) == 0) {
            j++;
            // debug_print("%s", str.s);
            ksprintf(&s, ">%d_%d%s\n", j, z.nsr, b->name);
            int j;
            for (j = 0; j < str.l; ++j) kputc("$ACGTN"[(int)str.s[j]], &s);
            kputc('\n', &s);
        }
    }
    */
    /*
    kstring_t s = {0,0,0};

    int i;
    for (i = 0; i < mg->v.n; ++i) {
        magv_t *v = &mg->v.a[i];
        ksprintf(&s, ">%d_%d%s\n", i, v->nsr, b->name);
        int j;
        for (j = 0; j < v->len; ++j) kputc("$ACGTN"[(int)v->seq[j]], &s);
        kputc('\n', &s);            
    }
    */
    mag_t *g = fml_fmi2mag(args.assem_opt, e);
    
    int n_utg;
    fml_utg_t *utg = fml_mag2utg(args.assem_opt, g);
    // fml_assemble(args.assem_opt, b->n, b->b, &n_utg);
    char *out = rend_utg(b->name, utg, n_utg);
    fml_utg_destroy(n_utg, utg);
    free(b->name);
    free(b);
    // read_block_destory(b);

    // free(str.s);
    // free(cov.s);
    // free(a.used); free(a.bend);
    return out;
    
  empty_block:
    read_block_destory(b);
    rld_destroy(e);
    return NULL;
}

static void write_out(void *s)
{
    if (s == NULL) args.filter_block++;
    else {
        args.assem_block++;
        fputs((char*)s, args.out);
        free(s);
    }
    // free(s);
}


int LFR_assem(int argc, char **argv)
{
    double t_real;
    t_real = realtime();

    if (parse_args(argc, argv)) return usage();

    hts_tpool *p = hts_tpool_init(args.n_thread);
    hts_tpool_process *q = hts_tpool_process_init(p, args.n_thread*2, 0);
    hts_tpool_result *r;
    
    for (;;) {
        //int nseq;
        struct read_block *b = read_block();
        if (b == NULL) break;
        //if (b->n == 1) continue;
        int block;
        do {
            block = hts_tpool_dispatch2(p, q, run_it, b, 1);
            if ((r = hts_tpool_next_result(q))) {
                void *s = hts_tpool_result_data(r);
                write_out(s);
                hts_tpool_delete_result(r, 0);
            }
        }
        while (block == -1);
    }

    hts_tpool_process_flush(q);
    
    while ((r = hts_tpool_next_result(q))) {
        void *s = hts_tpool_result_data(r);
        write_out(s);
        hts_tpool_delete_result(r, 0);
    }
    
    hts_tpool_process_destroy(q);
    hts_tpool_destroy(p);

    fprintf(stderr, "Filter block: %"PRIu64"\n", args.filter_block);
    fprintf(stderr, "Assembled block: %"PRIu64"\n", args.assem_block);
    memory_release();
    LOG_print("Real time: %.3f sec; CPU: %.3f sec", realtime() - t_real, cputime());
    return 0;
}
/*
int main(int argc, char **argv)
{
    return assem(argc, argv);
}
*/