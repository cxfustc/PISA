/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2021-03-12 13:04:59
  *Edit History: 
***********************************************************/

#include "str.h"

struct ss_s;
typedef struct ss_s ss_t;

ss_t * ss_init (int nth);

int ss_push (ss_t * engine, char * seq, int l);

str_t * ss_query (ss_t * engine, char * seq, int l, int dist, int * is_exact_match, int thid);

void ss_destroy (ss_t * engine);
