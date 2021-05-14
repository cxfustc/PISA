/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2021-03-12 14:44:01
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "str.h"
#include "hash.h"
#include "array.h"
#include "util.h"
#include "sim_search.h"

static int ss_kmer_size;

typedef struct {
  arr_t(i32) * idxs;
} sidx_t;

static void
sidx_init2 (sidx_t * si)
{
  si->idxs = arr_init (i32);
}

static void
sidx_clear (sidx_t * si)
{
  arr_clear (i32, si->idxs);
}

static void
sidx_copy (sidx_t * dst, sidx_t * src)
{
  arr_copy (i32, dst->idxs, src->idxs);
}

HASH_MAP_DEF (k2i, uint64_t, sidx_t);

struct ss_s {
  str_set_t * bcs;
  xh_map_t(u64_i32) * bc_hash;
  xh_map_t(k2i) * kmer_hash;

  sidx_t * si_buf;

  int nth;
  arr_t(i32) ** idx_sets;
};

static int
hasN (char * seq, int l)
{
  int i;

  for (i=0; i<l; ++i)
    if (base2int_tbl[seq[i]] == 4)
      return 1;

  return 0;
}

static void
build_kmers (xh_map_t(k2i) * hash, uint64_t cs, int32_t idx, sidx_t * si_buf)
{
  int offset;
  uint64_t kmer;
  uint64_t mask;
  sidx_t * si;

  offset = 3 * ss_kmer_size;
  mask = (1 << offset) - 1;
  for (;;) {
    if ((cs >> (offset-3)) == 0)
      break;

    kmer = cs & mask;
    cs = cs >> 3;
    if ((si = xh_map_search(k2i,hash,&kmer)) == NULL) {
      sidx_clear (si_buf);
      arr_add (i32, si_buf->idxs, idx);
      xh_map_add (k2i, hash, &kmer, si_buf);
    } else
      arr_add (i32, si->idxs, idx);
  }
}

ss_t *
ss_init (int nth)
{
  int i;
  ss_t * engine;

  engine = (ss_t *) ckmalloc (sizeof(ss_t));
  engine->bc_hash = xh_u64_i32_map_init ();
  engine->kmer_hash = xh_map_init (k2i, 256, 0.75,
      NULL, sidx_init2, NULL, sidx_copy, u32_hash_f, u32_equal_f);
  engine->bcs = str_set_init ();

  engine->si_buf = (sidx_t *) ckmalloc (sizeof(sidx_t));
  sidx_init2 (engine->si_buf);

  engine->nth = nth;
  engine->idx_sets = (arr_t(i32) **) ckmalloc (nth * sizeof(arr_t(i32) *));
  for (i=0; i<nth; ++i)
    engine->idx_sets[i] = arr_init (i32);

  ss_kmer_size = 5;

  return engine;
}

int
ss_push (ss_t * engine, char * seq, int l)
{
	int32_t idx;
	int32_t * ptr;
  uint64_t cs;

  if (hasN(seq,l))
    err_mesg ("Try to push sequence %s containing Ns.", seq);

  cs = enc64 (seq, l);

	if ((ptr = xh_map_search(u64_i32,engine->bc_hash,&cs)) == NULL) {
		idx = engine->bcs->n;
		xh_map_add (u64_i32, engine->bc_hash, &cs, &idx);
    str_set_add2 (engine->bcs, seq, l);
    build_kmers (engine->kmer_hash, cs, engine->bcs->n-1, engine->si_buf);
	}

	return 0;
}

static int
hamming_distance (str_t * ref, char * qry, int l)
{
  int i, d;

  for (i=d=0; i<l; ++i)
    if (ref->s[i] != qry[i])
      ++d;

  return d;
}

str_t *
ss_query (ss_t * engine, char * seq, int l, int dist, int * is_exact_match, int thid)
{
  int i;
  int d;
  int hit;
	int32_t * ptr;
  uint64_t cs;
  arr_t(i32) * set;
  sidx_t * si;

  cs = enc64 (seq, l);
  if ((ptr = xh_map_search(u64_i32,engine->bc_hash,&cs)) != NULL) {
    *is_exact_match = 1;
    return str_set_at (engine->bcs,*ptr);
  }

  set = engine->idx_sets[thid];
  arr_clear (i32, set);
  for (i=0; i<l-ss_kmer_size+1; ++i) {
    if (hasN(seq+i,ss_kmer_size)) {
      i += 1; // TODO: is this right?
      continue;
    }
    cs = enc64 (seq+i,ss_kmer_size);
    if ((si = xh_map_search(k2i,engine->kmer_hash,&cs)) == NULL)
      continue;
    arr_append (i32, set, si->idxs);
  }
  arr_uniqsort (i32, set);

  hit = -1;
  for (i=0; i<set->n; ++i) {
    d = hamming_distance (str_set_at(engine->bcs,set->arr[i]), seq, l);
    if (d <= dist) {
      if (hit != -1)
        return NULL;
      hit = set->arr[i];
    }
  }

  return str_set_at (engine->bcs,hit);
}

void
ss_destroy (ss_t * engine)
{

}

#if 0
int main()
{
    ss_t *S = ss_init(1);    
		printf ("ref\n");
    ss_push(S,"CTTCGATGGT", 10);
		printf ("\n");
    ss_push(S,"ACTTCTATGC", 10);
		printf ("\n");
    int exact;
		printf ("qry\n");
    str_t *s1 = ss_query(S, "CTTCTATGGT", 10, 1, &exact, 0);
		printf ("\n");
    str_t *s2 = ss_query(S, "ACTTCTATGA", 10, 2, &exact, 0);
		printf ("\n");
		if (s1)
    	printf("%s\n", s1->s);
		if (s2)
    	printf("%s\n", s2->s);
    return 0;
}
#endif
