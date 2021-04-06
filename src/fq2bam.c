// convert fastq file to raw bam file, no alignment info
#include "utils.h"
#include "number.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)
    
static struct args {
    const char *input_fname;
    const char *output_fname;
    gzFile rf;
    kseq_t *ks;
    int threads;
} args = {
    .input_fname = NULL,
    .output_fname = NULL,
    .rf = NULL,
    .ks = NULL,
    .threads =1,
};
static int usage()
{
    fprintf(stderr, "fq2bam -out out.bam in.fq\n");
    fprintf(stderr, "       -@  Threads to pack BAM file.\n");
    return 1;
}
static int parse_args(int argc, char **argv)
{
    if (argc == 1) return 1;
    const char *thrd = NULL;
    int i;
    const char **var = NULL;
    for (i = 1; i < argc; ) {
        const char *a = argv[i++];
        if (strcmp(a, "-h") == 0) return 1;
        if (strcmp(a, "-out") == 0 || strcmp(a, "-o") == 0) var = &args.output_fname;
        if (strcmp(a, "-@") == 0) var = &thrd;
        
        if (var != 0) {
            if (i == argc) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if (args.input_fname ==0) {
            args.input_fname = a;
            continue;
        }
        error("Unknown argument: %s, use -h see help infomation.", a);
    }

    args.rf = gzopen(args.input_fname, "r");
    if (args.rf == NULL) error("%s : %s.", args.input_fname, strerror(errno));

    args.ks = kseq_init(args.rf);

    if (thrd) {
        args.threads = str2int(thrd);
        if (args.threads < 1) args.threads = 1;
    }

    return 0;
}

int fq2bam(int argc, char **argv)
{
    if (parse_args(argc, argv)) return usage();

    htsFile *fp = hts_open(args.output_fname, "bw");
    if (fp == NULL) error("%s : %s.", args.output_fname, strerror(errno));
    const char *hdr_str = "@HD\tVN:1.4";
    int l_hdr_str = strlen(hdr_str);

    sam_hdr_t *hdr = sam_hdr_parse(l_hdr_str, hdr_str);

    if (sam_hdr_write(fp, hdr) != 0) error("Failed to write hdr.");

    hts_set_threads(fp, args.threads);
    
    kstring_t str = {0,0,0};
    kstring_t temp = {0,0,0};
    while(1) {
        str.l = 0;
        temp.l = 0;
        if (kseq_read(args.ks) <0) break;
        kstring_t *name = &args.ks->name;
        kstring_t *seq  = &args.ks->seq;
        kstring_t *qual = &args.ks->qual;

        int i = 0;
        for (i = 0; i < name->l && name->s[i]!='|'; ++i);
        if (i > 0 && i < name->l -5) {
            char *p = name->s+i;
            char *r = 0;
            char *e = name->s+name->l;
            *e = '\0';
            for (;p != e;) {
                if (*p == '|' && *(p+1) == '|' && *(p+2) == '|') {
                    *p = 0;
                    if (r != 0) {
                        kputc('\t', &temp);
                        kputs(r, &temp);
                    }
                    p+=3;
                    r=p;
                }
                p++;
            }
            if (r) {
                kputc('\t', &temp);
                kputs(r, &temp);
            }
        }

        kputs(name->s, &str);

        kputs("\t4\t*\t0\t0\t*\t*\t0\t0\t", &str);
        kputs(seq->s, &str);
        kputc('\t', &str);
        if (qual->l == 0) kputc('*', &str);
        else kputs(qual->s, &str);

        kputs(temp.s, &str);
        
        bam1_t *b = bam_init1();
        if (sam_parse1(&str, hdr, b)) error("Failed to parse.");
        if (sam_write1(fp, hdr, b) == -1) error("Failed to write.");
        bam_destroy1(b);
    }

    if (temp.l) free(temp.s);
    if (str.l) free(str.s);

    kseq_destroy(args.ks);
    gzclose(args.rf);

    hts_close(fp);
    return 0;
}
