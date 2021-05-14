/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2021-03-18 14:55:34
  *Edit History: 
***********************************************************/

#ifndef XDK_KSEQ_H
#define XDK_KSEQ_H

#include <kseq.h>

#include "rseq.h"
#include "util.h"

KSEQ_INIT(gzFile, err_gzread)

struct pefq_kseq_reader_s;
typedef struct pefq_kseq_reader_s pefq_kseq_reader_t;

pefq_kseq_reader_t * pefq_kseq_reader_init (const char * fq1_file, const char * fq2_file);
void pefq_kseq_reader_destroy (pefq_kseq_reader_t * reader);
int pefq_kseq_reader_load_block (pefq_kseq_reader_t * reader, pefq_blk_t * blk);

#endif
