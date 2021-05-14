/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2018-12-10 17:34:56
  *Edit History: 
***********************************************************/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>

#include "str.h"
#include "rseq.h"
#include "util.h"

#define INIT_SIZE 128
#define INIT_SIZE_BT_WIDTH 7

/*****************************************************************************/
/*------------------------- read seq definition -----------------------------*/
/*****************************************************************************/

rseq_t *
rseq_init (void)
{
	rseq_t * r;

	r = (rseq_t *) ckalloc (1, sizeof(rseq_t));
	rseq_init2 (r);

	return r;
}

void
rseq_init2 (rseq_t * r)
{
	r->l = 0;
	r->m = INIT_SIZE;
	r->b = (char *) ckmalloc (INIT_SIZE);
	r->q = (char *) ckmalloc (INIT_SIZE);
}

void
rseq_resize (rseq_t * r, int32_t new_size)
{
	if (new_size < r->m)
		return;

	r->m = ((new_size>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
	free (r->b);
	free (r->q);
	r->b = (char *) ckmalloc (r->m);
	r->q = (char *) ckmalloc (r->m);
}

void
rseq_clear (rseq_t * r)
{
	r->l = 0;
}

void
rseq_free (rseq_t * r)
{
	rseq_free2 (r);
	free (r);
}

void
rseq_free2 (rseq_t * r)
{
	free (r->b);
	free (r->q);
}

void
rseq_assign (rseq_t * r, const char * name, const char * anno,
                const char * bases, const char * quals)
{
  int32_t l_seq;

  l_seq = strlen (bases);
  if (l_seq != strlen(quals))
    err_mesg ("bases and quals are not same length!");

  rseq_resize (r, l_seq);
  memcpy (r->b, bases, l_seq);
  memcpy (r->q, quals, l_seq);
  r->l = l_seq;
}

void
rseq_reverse_complement (rseq_t * r)
{
  char tmp;
  int32_t i, j;
  int32_t mid;

  mid = r->l >> 1;
  for (i=0,j=r->l-1; i<mid; ++i,--j) {
    tmp = r->b[i];
    r->b[i] = base_rc_tbl[r->b[j]];
    r->b[j] = base_rc_tbl[tmp];

    tmp = r->q[i];
    r->q[i] = r->q[j];
    r->q[j] = tmp;
  }

  if (r->l & 0x1)
    r->b[mid] = base_rc_tbl[r->b[mid]];
}

/*****************************************************************************/
/*---------------------- pair end read seq definition -----------------------*/
/*****************************************************************************/

pe_seq_t *
pe_seq_init (void)
{
	pe_seq_t * p;

	p = (pe_seq_t *) ckmalloc (sizeof(pe_seq_t));
	pe_seq_init2 (p);

	return p;
}

void
pe_seq_init2 (pe_seq_t * p)
{
  p->n = str_init ();
	p->l[0] = p->l[1] = 0;
	p->m[0] = p->m[1] = INIT_SIZE;
  p->b[0] = ckmalloc (p->m[0]);
  p->b[1] = ckmalloc (p->m[0]);
  p->q[0] = ckmalloc (p->m[0]);
  p->q[1] = ckmalloc (p->m[0]);
}

void
pe_seq_resize (pe_seq_t * p, int32_t new_size1, int32_t new_size2)
{
  if (new_size1 >= p->m[0]) {
	  p->m[0] = ((new_size1>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
	  free (p->b[0]);
	  free (p->q[0]);
	  p->b[0] = ckmalloc (p->m[0]);
	  p->q[0] = ckmalloc (p->m[0]);
  }

  if (new_size2 >= p->m[1]) {
	  p->m[1] = ((new_size2>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
	  free (p->b[1]);
	  free (p->q[1]);
	  p->b[1] = ckmalloc (p->m[1]);
	  p->q[1] = ckmalloc (p->m[1]);
  }
}

void
pe_seq_clear (pe_seq_t * p)
{
  p->l[0] = p->l[1] = 0;
  str_clear (p->n);
}

void
pe_seq_free (pe_seq_t * p)
{
  pe_seq_free2 (p);
  free (p);
}

void
pe_seq_free2 (pe_seq_t * p)
{
  free (p->b[0]);
  free (p->b[1]);
  free (p->q[0]);
  free (p->q[1]);
  str_free (p->n);
}

void
pe_seq_assign (pe_seq_t * p, int end, const char * name, const char * anno,
    const char * bases, const char * quals)
{

}

void
pe_seq_pool_dump (FILE ** out, mp_t(pe) * pes)
{
  int64_t i, j;
  pe_seq_t * r;

  for (i=0; i<mp_cnt(pes); ++i) {
    r = mp_at (pe, pes, i);
    for (j=0; j<2; ++j) {
      fprintf (out[j], "%s/1\n", r->n->s);
      fprintf (out[j], "%s\n", r->b[j]);
      fprintf (out[j], "+\n");
      fprintf (out[j], "%s\n", r->q[j]);
    }
  }
}

void
pe_seq_pool_dump_lfq (FILE * out, mp_t(pe) * pes)
{
  int64_t i;
  pe_seq_t * r;

  for (i=0; i<mp_cnt(pes); ++i) {
    r = mp_at (pe, pes, i);
    fprintf (out, "%s %s\n", r->b[0], r->b[1]);
  }
}

/*****************************************************************************/
/*-------------------------- block of read seqs -----------------------------*/
/*****************************************************************************/

#define RSEQ_BLK_SIZE 4096

void
rseq_blk_init2 (rseq_blk_t * blk)
{
  int32_t i;

  blk->n = 0;
  blk->m = RSEQ_BLK_SIZE;
  blk->seqs = (rseq_t *) ckmalloc (blk->m * sizeof(rseq_t));
  for (i=0; i<RSEQ_BLK_SIZE; ++i)
    rseq_init2 (blk->seqs+i);
}

void
rseq_blk_clear (rseq_blk_t * blk)
{
  blk->n = 0;
}

void
rseq_blk_free2 (rseq_blk_t * blk)
{
  int i;

  for (i=0; i<blk->m; ++i)
    rseq_free (blk->seqs+i);
  free (blk->seqs);
}

/*****************************************************************************/
/*--------------------------- single fastq reader ---------------------------*/
/*****************************************************************************/

/*
static void *
sefq_reader_core (void * data)
{
  char ch;
  char * buf;
  int nth;
  int32_t os;
  int64_t n_line;
  gzFile fp;
  rseq_t * r;
  rseq_blk_t ** blks;
  sefq_reader_t * reader;

  reader = (sefq_reader_t *) data;
  nth = reader->nth;
  blks = reader->blks;

  n_line = os = 0;
  if ((fp = gzopen(reader->file,"r")) == NULL)
    err_mesg ("fail to open file: '%s'!", reader->file);

  r = blks[0][0].seqs + 0;

  while ((ch = gzgetc(fp)) != -1) {
    if (ch == '\n') {
      if ((n_line & 0x3) == 0) {
        buf
      }
    }
  }
}

sefq_reader_t *
sefq_reader_init (const char * fq_file, int nth)
{
  int32_t i;
  sefq_reader_t * reader;

  reader = (sefq_reader_t *) ckmalloc (sefq_reader_t);
  reader->nth = nth;
  reader->file = strdup (fq_file);
  reader->blks = (rseq_blk_t **) ckmalloc (nth * sizeof(sefq_reader_t *));
  for (i=0; i<nth; ++i) {
    reader->blks[i] = (rseq_blk_t *) ckmalloc (2 * sizeof(sefq_reader_t));
    rseq_blk_init2 (reader->blks[i]);
    rseq_blk_init2 (reader->blks[i]+1);
  }

  ckpthread_create (&reader->pid, NULL, sefq_reader_core, (void*)(reader));

  return reader
}

void
sefq_reader_destroy (sefq_reader_t * reader)
{
  int32_t i;

  ckpthread_join (reader->pid);

  for (i=0; i<reader->nth; ++i) {
    rseq_blk_free2 (reader->blks[i]);
    rseq_blk_free2 (reader->blks[i]+1);
    free (reader->blks[i]);
  }
  free (reader->blks);
  free (reader);
}
*/

mp_t(rs) *
sefq_load (const char * fq_file)
{
  char ch;
  int32_t i;
  int64_t n_line;
  time_t time_beg;
  str_t * b;
  str_t * q;
  gzFile fp;
  rseq_t * r;
  mp_t(rs) * set;

  time (&time_beg);

  if ((fp = ckgzopen(fq_file,"r")) == 0)
    err_mesg ("fail to open file: '%s'!", fq_file);

  set = mp_init (rs, NULL);
  b = str_init ();
  q = str_init ();
  n_line = 0;

  while ((ch = gzgetc(fp)) != -1) {
    if (ch == '\n') {
      if ((n_line & 0x3) == 3) {
        if (b->l != q->l)
          err_mesg ("l_base != l_qual");
        r = mp_alloc (rs, set);
        r->l = b->l;
        r->m = ((r->l>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
        r->b = (char *) ckmalloc (r->m);
        r->q = (char *) ckmalloc (r->m);
        memcpy (r->b, b->s, b->l);
        memcpy (r->q, q->s, b->l);

        str_clear (b);
        str_clear (q);
      }
      ++n_line;
      continue;
    }

    if ((n_line & 0x3) == 1)
      str_add (b, ch);
    else if ((n_line & 0x3) == 3)
      str_add (q, ch);
  }

  gzclose (fp);

  str_free (b);
  //str_free (q);

  //for (i=0; i<mp_cnt(set); ++i) {
    //r = mp_at (rs, set, i);
    //printf ("> read length %d\n", r->l);
    //printf ("%s\n", r->b);
    //printf ("%s\n", r->q);
  //}

  printf ("Total loading %ld reads, cost %lds\n", mp_cnt(set), time(NULL)-time_beg);

  return set;
}

/*****************************************************************************/
/*--------------------------- paired fastq reader ---------------------------*/
/*****************************************************************************/

struct pefq_reader_s {
  gzFile fq[2];
  str_t * base;
  str_t * qual;
  str_t * name;
};

pefq_blk_t *
pefq_blk_init (int32_t max)
{
  int32_t i;
  pefq_blk_t * blk;

  blk = (pefq_blk_t *) ckalloc (1, sizeof(pefq_blk_t));
  blk->n = 0;
  blk->m = max;
  blk->seqs = (pe_seq_t *) ckalloc (max, sizeof(pe_seq_t));
  for (i=0; i<max; ++i)
    pe_seq_init2 (blk->seqs+i);

  return blk;
}

void
pefq_blk_free (pefq_blk_t * blk)
{
  int32_t i;

  for (i=0; i<blk->m; ++i)
    pe_seq_free2 (blk->seqs+i);
  free (blk->seqs);
  free (blk);
}

pefq_reader_t *
pefq_reader_init (const char * fq1_file, const char * fq2_file)
{
  pefq_reader_t * reader;

  reader = (pefq_reader_t *) ckalloc (1, sizeof(pefq_reader_t));
  reader->fq[0] = ckgzopen (fq1_file, "r");
  reader->fq[1] = ckgzopen (fq2_file, "r");
  reader->base = str_init ();
  reader->qual = str_init ();
  reader->name = str_init ();

  return reader;
}

void
pefq_reader_destroy (pefq_reader_t * reader)
{
  gzclose (reader->fq[0]);
  gzclose (reader->fq[1]);
  str_free (reader->base);
  str_free (reader->qual);
  str_free (reader->name);
  free (reader);
}

int
pefq_reader_load_block (pefq_reader_t * reader, pefq_blk_t * blk)
{
  char ch;
	time_t beg;
  int32_t n_lines;
  int32_t n_reads;
  gzFile fp;
  str_t * b;
  str_t * q;
  str_t * n;
  pe_seq_t * pe;

	time (&beg);

  b = reader->base;
  q = reader->qual;
  n = reader->name;

  blk->n = 0;
  n_lines =  0;
  str_clear (b);
  str_clear (q);
  str_clear (n);
  fp = reader->fq[0];
  while ((ch = gzgetc(fp)) != -1) {
    if (ch == '\n') {
      if ((n_lines & 0x3) == 3) {
        if (b->l != q->l)
          err_mesg ("pe reads, l_base != l_qual");

        pe = blk->seqs + blk->n++;
        pe->l[0] = b->l;
        pe_seq_resize (pe, b->l, -1);
        memcpy (pe->b[0], b->s, b->l);
        memcpy (pe->q[0], q->s, q->l);
        pe->b[0][b->l] = '\0';
        pe->q[0][q->l] = '\0';
				if (n->s[n->l-2]=='/' && n->s[n->l-1]=='1')
					n->s[n->l-2]=='\0', n->l-=2;
        str_copy (pe->n, n);
        pe->flag = 0;

        str_clear (b);
        str_clear (q);
        str_clear (n);
      }
      ++n_lines;
      if (blk->n == blk->m)
        break;
      continue;
    }

    if ((n_lines & 0x3) == 0)
      str_add (n, ch);
    else if ((n_lines & 0x3) == 1)
      str_add (b, ch);
    else if ((n_lines & 0x3) == 3)
      str_add (q, ch);
  }
  n_reads = blk->n;

  blk->n = 0;
  n_lines = 0;
  str_clear (b);
  str_clear (q);
  fp = reader->fq[1];
  while ((ch = gzgetc(fp)) != -1) {
    if (ch == '\n') {
      if ((n_lines & 0x3) == 3) {
        if (b->l != q->l)
          err_mesg ("pe reads, l_base != l_qual");

        pe = blk->seqs + blk->n++;
        pe->l[1] = b->l;
        pe_seq_resize (pe, -1, b->l);
        memcpy (pe->b[1], b->s, b->l);
        memcpy (pe->q[1], q->s, q->l);
        pe->b[1][b->l] = '\0';
        pe->q[1][q->l] = '\0';

        str_clear (b);
        str_clear (q);
      }
      ++n_lines;
      if (blk->n == blk->m)
        break;
      continue;
    }

    if ((n_lines & 0x3) == 1)
      str_add (b, ch);
    else if ((n_lines & 0x3) == 3)
      str_add (q, ch);
  }

  assert (n_reads == blk->n);

	printf ("Fastq block loading costs: %lds\n", time(NULL)-beg);

  if (blk->n > 0)
    return 0;
  else
    return -1;
}

mp_t(pe) *
pefq_load (const char * fq1_file, const char * fq2_file)
{
  char ch;
  int64_t rd_idx;
  int64_t rd_cnt;
  int64_t n_line;
  time_t time_beg;
  str_t * b;
  str_t * q;
  str_t * n;
  gzFile fp;
  pe_seq_t * pe;
  mp_t(pe) * set;

  time (&time_beg);

  set = mp_init (pe, NULL);
  b = str_init ();
  q = str_init ();
  n = str_init ();
  n_line = 0;

  if ((fp = ckgzopen(fq1_file,"r")) == 0)
    err_mesg ("fail to open file: '%s'!", fq1_file);

  while ((ch = gzgetc(fp)) != -1) {
    if (ch == '\n') {
      if ((n_line & 0x3) == 3) {
        if (b->l != q->l)
          err_mesg ("lfr reads, l_base != l_qual");

        pe = mp_alloc (pe, set);
        pe->l[0] = b->l;
        pe->m[0] = ((pe->l[0]>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
        pe->b[0] = (char *) ckmalloc (pe->m[0]);
        pe->q[0] = (char *) ckmalloc (pe->m[0]);
        memcpy (pe->b[0], b->s, b->l);
        memcpy (pe->q[0], q->s, q->l);
        pe->b[0][b->l] = '\0';
        pe->q[0][q->l] = '\0';
        pe->n = str_init ();
        str_copy (pe->n, n);

        str_clear (b);
        str_clear (q);
        str_clear (n);
      }
      ++n_line;
      continue;
    }

    if ((n_line & 0x3) == 0)
      str_add (n, ch);
    else if ((n_line & 0x3) == 1)
      str_add (b, ch);
    else if ((n_line & 0x3) == 3)
      str_add (q, ch);
  }
  gzclose (fp);

  if ((fp = ckgzopen(fq2_file,"r")) == 0)
    err_mesg ("fail to open file: '%s'!", fq2_file);

  str_clear (b);
  str_clear (q);
  rd_idx = 0;
  rd_cnt = mp_cnt (set);
  while ((ch = gzgetc(fp)) != -1) {
    if (ch == '\n') {
      if ((n_line & 0x3) == 3) {
        if (b->l != q->l)
          err_mesg ("lfr reads, l_base != l_qual");

        if (rd_idx >= rd_cnt)
          err_mesg ("there are diffent number of reads in '%s' and '%s'!", fq1_file, fq2_file);

        pe = mp_at (pe, set, rd_idx++);
        pe->l[1] = b->l;
        pe->m[1] = ((pe->l[1]>>INIT_SIZE_BT_WIDTH)+1) << INIT_SIZE_BT_WIDTH;
        pe->b[1] = (char *) ckmalloc (pe->m[1]);
        pe->q[1] = (char *) ckmalloc (pe->m[1]);
        memcpy (pe->b[1], b->s, b->l);
        memcpy (pe->q[1], q->s, q->l);
        pe->b[1][b->l] = '\0';
        pe->q[1][q->l] = '\0';

        str_clear (b);
        str_clear (q);
      }
      ++n_line;
      continue;
    }

    if ((n_line & 0x3) == 1)
      str_add (b, ch);
    else if ((n_line & 0x3) == 3)
      str_add (q, ch);
  }
  gzclose (fp);

  str_free (n);
  str_free (b);
  str_free (q);

  if (rd_idx != rd_cnt)
    err_mesg ("there are diffent number of reads in '%s' and '%s'!", fq1_file, fq2_file);

  printf ("'%s' and '%s':\n", fq1_file, fq2_file);
  printf ("Total loading %ld pair end reads, cost %lds\n",
      rd_cnt, time(NULL)-time_beg);

  return set;
}

void
pefq_load1pair (mp_t(pe) * set, const char * fq1_file, const char * fq2_file)
{

  char ch;
  int64_t new_rd_idx;
  int64_t new_rd_cnt;
  int64_t old_rd_cnt;
  int64_t n_line;
  time_t time_beg;
  str_t * b;
  str_t * q;
  str_t * n;
  gzFile fp;
  pe_seq_t * pe;

  time (&time_beg);

  b = str_init ();
  q = str_init ();
  n = str_init ();
  n_line = 0;

  old_rd_cnt = mp_cnt (set);
  if ((fp = ckgzopen(fq1_file,"r")) == 0)
    err_mesg ("fail to open file: '%s'!", fq1_file);

  while ((ch = gzgetc(fp)) != -1) {
    if (ch == '\n') {
      if ((n_line & 0x3) == 3) {
        if (b->l != q->l)
          err_mesg ("lfr reads, l_base != l_qual");

        pe = mp_alloc (pe, set);
        pe->l[0] = b->l;
        pe_seq_resize (pe, b->l, -1);
        memcpy (pe->b[0], b->s, b->l);
        memcpy (pe->q[0], q->s, q->l);
        pe->b[0][b->l] = '\0';
        pe->q[0][q->l] = '\0';
        str_copy (pe->n, n);

        str_clear (b);
        str_clear (q);
        str_clear (n);
      }
      ++n_line;
      continue;
    }

    if ((n_line & 0x3) == 0)
      str_add (n, ch);
    else if ((n_line & 0x3) == 1)
      str_add (b, ch);
    else if ((n_line & 0x3) == 3)
      str_add (q, ch);
  }
  gzclose (fp);

  if ((fp = ckgzopen(fq2_file,"r")) == 0)
    err_mesg ("fail to open file: '%s'!", fq2_file);

  str_clear (b);
  str_clear (q);
  new_rd_idx = old_rd_cnt;
  new_rd_cnt = mp_cnt (set);
  while ((ch = gzgetc(fp)) != -1) {
    if (ch == '\n') {
      if ((n_line & 0x3) == 3) {
        if (b->l != q->l)
          err_mesg ("lfr reads, l_base != l_qual");

        if (new_rd_idx >= new_rd_cnt)
          err_mesg ("there are diffent number of reads in '%s' and '%s'!", fq1_file, fq2_file);

        pe = mp_at (pe, set, new_rd_idx++);
        pe->l[1] = b->l;
        pe_seq_resize (pe, -1, b->l);
        memcpy (pe->b[1], b->s, b->l);
        memcpy (pe->q[1], q->s, q->l);
        pe->b[1][b->l] = '\0';
        pe->q[1][q->l] = '\0';

        str_clear (b);
        str_clear (q);
      }
      ++n_line;
      continue;
    }

    if ((n_line & 0x3) == 1)
      str_add (b, ch);
    else if ((n_line & 0x3) == 3)
      str_add (q, ch);
  }
  gzclose (fp);

  str_free (n);
  str_free (b);
  str_free (q);

  if (new_rd_idx != new_rd_cnt)
    err_mesg ("there are diffent number of reads in '%s' and '%s'!", fq1_file, fq2_file);

  //printf ("'%s' and '%s':\n", fq1_file, fq2_file);
  //printf ("Total loading %ld pair end reads, cost %lds\n",
  //    new_rd_cnt-old_rd_cnt, time(NULL)-time_beg);
}
