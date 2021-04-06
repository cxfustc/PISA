#ifndef GTF_H
#define GTF_H

#include <stdlib.h>
#include "dict.h"
#include "region_index.h"
#define EXONIC    0
#define INTRONIC  1

enum feature_type {
    feature_unknow =-1,
    feature_gene,
    feature_transcript,
    feature_lncRNA,
    feature_CDS,
    feature_start_codon,
    feature_stop_codon,
    feature_5UTR,
    feature_3UTR,
    feature_inter,
    feature_inter_CNS,
    feature_intron_CNS,
    feature_exon,
    feature_5UTR_alias,
    feature_3UTR_alias,
    feature_Selenocysteine,
};

static const char *feature_type_names[] = {
    // The following feature types are required: "gene", "transcript"
    "gene",
    "transcript",
    "lnc_RNA",
    // The features "CDS", "start_codon", "stop_codon", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are optional.
    "CDS",
    "start_codon",
    "stop_codon",
    "5UTR",
    "3UTR",
    "inter",
    "inter_CNS",
    "intron_CNS",
    "exon",
    "five_prime_utr",
    "three_prime_utr",
    "Selenocysteine"
    // All other features will be ignored. The types must have the correct capitalization shown here.
};
struct gtf {
    int seqname;
    int source;
    enum feature_type type;
    int start;
    int end;
    int strand; // 0 on forward, 1 on reverse
    int gene_id;
    int gene_name;
    int transcript_id;

    struct dict *attr; // attributions
    struct dict *query; // used to fast access gtf
    int n_gtf, m_gtf;
    struct gtf **gtf;
};

struct _ctg_idx;

struct gtf_ctg {
    struct dict *gene_idx;
    struct region_index *idx;    
    int n_gtf, m_gtf;
    struct gtf **gtf; 
};
struct gtf_spec {
    struct dict *name; // contig names
    struct dict *gene_name;
    struct dict *gene_id;
    struct dict *transcript_id;
    struct dict *sources; //
    struct dict *attrs; // attributes
    struct dict *features;
};

const char *get_feature_name(enum feature_type type);

struct gtf_spec *gtf_read(const char *fname, int filter);
struct gtf_spec *gtf_read_lite(const char *fname); // only read necessary info
struct region_itr *gtf_query(struct gtf_spec const *G, char *name, int start, int end);
void gtf_destroy(struct gtf_spec *G);

#endif
