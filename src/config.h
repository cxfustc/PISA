/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2021-03-11 14:47:18
  *Edit History: 
***********************************************************/

#ifndef XDK_CONFIG_H
#define XDK_CONFIG_H

#include "sim_search.h"

struct config_s;
typedef struct config_s config_t;

struct bcode_reg {
    int rd; // read 1 or 2
    int start;
    int end;
    int dist;
    ss_t *wl;
    char **white_list; // temp allocated, will be free after initization
    int len;
    int n_wl;
};

struct config_s {
    // consistant with white list if set
    char *cell_barcode_tag; 
    char *sample_barcode_tag;

    char *raw_cell_barcode_tag;
    char *raw_cell_barcode_qual_tag;
    int n_cell_barcode;
    struct bcode_reg *cell_barcodes;

    char *raw_sample_barcode_tag;
    char *raw_sample_barcode_qual_tag;
    int n_sample_barcode; // usually == 1, sometimes == 2,
    struct bcode_reg *sample_barcodes;

    char *umi_tag;
    char *umi_qual_tag;
    // UMI, usually random generated and located at one region
    struct bcode_reg *UMI;

    // clean read sequence
    struct bcode_reg *read_1; 
    struct bcode_reg *read_2;  // for single ends, read2 == NULL
};

config_t * config_init2 (const char *fn, int nth);
void config_destroy2 (config_t * config);

#endif
