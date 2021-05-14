/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2021-03-11 12:47:13
  *Edit History: 
***********************************************************/

#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include "str.h"
#include "xth.h"
#include "hash.h"
#include "rseq.h"
#include "utils.h"
#include "config.h"
#include "str_hash.h"
#include "math_basic.h"

#define SIG_WAIT 0
#define SIG_WORK 1
#define SIG_STOP 2

#define BATCH_RD_CNT 1000000

#define FQ_FLAG_READ_LOWQ  8
#define PE_FLAG_SBC_FAIL   16
#define PE_FLAG_CBC_FAIL   32
#define PE_CBC_EXACT_MATCH 64

#define OPT_FQ1_OUT_GZ    1
#define OPT_FQ2_OUT_GZ    2

typedef struct {
  char * out_file[2];
  char * cfg_file;
  char * cbdis_file;
  char * report_file;
  int16_t min_qual;
  int16_t dropN;
	int32_t bgiseq_filter;

  char ** in_file[2];
  int16_t nth;
  int16_t n;

  uint32_t flag;
  FILE * fp[2];
  gzFile gfp[2];
} opt_t;

typedef struct {
  uint64_t raw_reads;
  uint64_t reads_pass_qc;
  uint64_t bc_exact_match;
  uint64_t filter_by_bc;
  uint64_t filter_by_lowqual;
  uint64_t filter_by_sample;

  uint64_t * bases_cell_barcode;
  uint64_t * q30_bases_cell_barcode;
  uint64_t * bases_sample_barcode;
  uint64_t * q30_bases_sample_barcode;
  uint64_t * bases_umi;
  uint64_t * q30_bases_umi;
  uint64_t * bases_read;
  uint64_t * q30_bases_read;
} qual_stat_t;

typedef struct {
  uint64_t bases_cell_barcode;
  uint64_t q30_bases_cell_barcode;
  uint64_t bases_sample_barcode;
  uint64_t q30_bases_sample_barcode;
  uint64_t bases_umi;
  uint64_t q30_bases_umi;
  uint64_t bases_read;
  uint64_t q30_bases_read;
} read_stat_t;

static opt_t * opt;
static config_t * cfg;
static pefq_blk_t * blk;
static pefq_reader_t * reader;
static int32_t cnter;
//static str_t * cbc_buf;
//static str_t * read_buf;
static str_hash_t ** bc_hashs;
static qual_stat_t st;

static int
usage (void)
{
  fprintf (stderr, "\n");
  fprintf (stderr, "Parse cell barcode and UMI string from raw FASTQ.\n");
  fprintf (stderr, "scRNA_v2_parse [options] lane1_1.fq.gz,lane2_1.fq.gz lane1_2.fq.gz,lane2_2.fq.gz\n");
  fprintf (stderr, "version: v 1.0.0\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Options :\n");
  fprintf (stderr, " -1       [fastq]   Read 1 output.\n");
  fprintf (stderr, " -2       [fastq]   Read 2 output.\n");
  fprintf (stderr, " -f                 Process bgiseq filter.\n");
  fprintf (stderr, " -t       [INT]     number of threads.\n");
  fprintf (stderr, " -config  [json]    Configure file in JSON format. Required.\n");
  //fprintf (stderr, " -run     [string]  Run code, used for different library.\n");
  fprintf (stderr, " -cbdis   [file]    Read count per cell barcode.\n\n");
  fprintf (stderr, " -q       [INT]     Drop reads if average sequencing quality below this value.\n");
  fprintf (stderr, " -dropN             Drop reads if N base in sequence or barcode.\n");
  fprintf (stderr, " -report  [csv]     summary report\n");
  fprintf (stderr, "\n");

  return 1;
}

static void
dump_options (opt_t * opt)
{
  int i;

  printf ("Input fastqs:\n");
  for (i=0; i<opt->n; ++i) {
    printf ("  ------ lane %d ------\n", i+1);
    printf ("  %s\n", opt->in_file[0][i]);
    printf ("  %s\n", opt->in_file[1][i]);
  }
  printf ("\n");

  printf ("Input config JSON file:\n");
  printf ("  %s\n", opt->cfg_file);
  printf ("\n");

  printf ("Output read count file:\n");
  printf ("  %s\n", opt->cbdis_file);
  printf ("\n");

  printf ("Output report file:\n");
  printf ("  %s\n", opt->report_file);
  printf ("\n");

  printf ("Minimum base quality: %d\n", opt->min_qual);
  printf ("DropN: %s\n", opt->dropN ? "True" : "False");
  printf ("\n");
}

static int
fq_file_count (char * file_list)
{
  char * ch;
  int n;

  for (n=1,ch=file_list; *ch; ++ch)
    if (*ch == ',')
      ++n;

  return n;
}

static void
save_file_names (char ** fnames, int n, char * list)
{
  int i = 0;
  char * ch;

  ch = strtok (list, ",");
  fnames[i++] = strdup (ch);

  while ((ch = strtok(NULL,",")) != NULL) {
    assert (i < n);
    fnames[i++] = strdup (ch);
  }

  assert (i == n);
}

static void
open_out_file (opt_t * opt, int idx)
{
  if (opt->out_file[idx] == NULL)
    return;

  if (str_cmp_tail(opt->out_file[idx],".fastq",6) == 0
      || str_cmp_tail(opt->out_file[idx],".fq",3) == 0) {
    opt->fp[idx] = ckopen (opt->out_file[idx], "w");
  } else if (str_cmp_tail(opt->out_file[idx],".fastq.gz",9) == 0
      || str_cmp_tail(opt->out_file[idx],".fq.gz",6) == 0) {
    opt->flag |= 1 << idx;
    opt->gfp[idx] = ckgzopen (opt->out_file[idx], "w");
  }
}

static opt_t *
parse_options (int argc, char * argv[])
{
  if (argc < 3) {
    usage ();
    abort ();
  }

  int i = 0;
  int n1, n2;
  char * key;
  opt_t * opt;

  opt = (opt_t *) ckalloc (1, sizeof(opt_t));
  opt->nth = 1;
  opt->flag = 0;

  while (i < argc-2) {
    key = argv[i++];

    if (strcmp(key,"-h")==0 || strcmp(key,"--help")==0 || strcmp(key,"-help")==0) {
      usage ();
      abort ();
    } else if (strcmp(key, "-1") == 0)
      opt->out_file[0] = strdup (argv[i++]);
    else if (strcmp(key, "-2") == 0)
      opt->out_file[1] = strdup (argv[i++]);
    else if (strcmp(key, "-config") == 0)
      opt->cfg_file = strdup (argv[i++]);
    else if (strcmp(key, "-cbdis") == 0)
      opt->cbdis_file = strdup (argv[i++]);
    else if (strcmp(key, "-q") == 0)
      opt->min_qual = atoi (argv[i++]);
    else if (strcmp(key, "-dropN") == 0)
      opt->dropN = 1;
    else if (strcmp(key, "-report") == 0)
      opt->report_file = strdup (argv[i++]);
    else if (strcmp(key, "-t") == 0)
      opt->nth = atoi (argv[i++]);
		else if (strcmp(key, "-f") == 0)
			opt->bgiseq_filter = 1;
  }

  n1 = fq_file_count (argv[argc-2]);
  n2 = fq_file_count (argv[argc-1]);
  assert (n1==n2 && n1>0);

  opt->n = n1;
  opt->in_file[0] = (char **) ckalloc (n1, sizeof(char *));
  opt->in_file[1] = (char **) ckalloc (n2, sizeof(char *));
  save_file_names (opt->in_file[0], n1, argv[argc-2]);
  save_file_names (opt->in_file[1], n2, argv[argc-1]);

  open_out_file (opt, 0);
  open_out_file (opt, 1);

  dump_options (opt);

  return opt;
}

static void
read_name_update (str_t * name, char * tag, char * s, int l)
{
  str_append (name, "|||", 3);
  str_append (name, tag, strlen(tag));
  str_append (name, ":Z:", 3);
  str_append (name, s, l);
}

static void
umi_process (pe_seq_t * seq, int thid, read_stat_t * rs)
{
  char * b;
  char * q;
  int i, l;
  int ridx;
  int start;
  int end;
  char * tag;
  char * qual_tag;

  ridx = cfg->UMI->rd-1;
  start = cfg->UMI->start;
  end = cfg->UMI->end;
  tag = cfg->umi_tag;
  qual_tag = cfg->umi_qual_tag;

  b = seq->b[ridx] + start - 1;
  q = seq->q[ridx] + start - 1;
  l = end - start + 1;

  rs->bases_umi += l;
  for (i=0; i<l; ++i)
    if (q[i] >= 63) // q30: 33 + 30
      ++(rs->q30_bases_umi);

  if (opt->dropN) {
    for (i=0; i<l; ++i)
      if (b[i] == 'N') {
        seq->flag = FQ_FLAG_READ_LOWQ;
        break;
      }
  }

  if (tag)
    read_name_update (seq->n, tag, b, l);

  if (qual_tag)
    read_name_update (seq->n, qual_tag, q, l);
}

static void
barcode_tag_create (char * b, char * q, int l, str_t * wl, str_t * buf)
{
  if (cfg->raw_cell_barcode_tag)
    str_append (buf+1, b, l);

  if (cfg->raw_cell_barcode_qual_tag)
    str_append (buf+2, q, l);

  if (cfg->cell_barcode_tag) {
    if (wl != NULL)
      b = wl->s;
    str_append (buf, b, l);
  }
}

// -1: filtered
//  0: others
static int
process1barcode (pe_seq_t * seq, struct bcode_reg * br, int thid, str_t * buf, int * em, read_stat_t * rs)
{
  int i;
  int l;
  char * b;
  char * q;
  str_t * s;

  b = seq->b[br->rd-1] + br->start - 1;
  q = seq->q[br->rd-1] + br->start - 1;
  l = br->end - br->start + 1;
  assert (l == br->len);

  rs->bases_cell_barcode += l;
  for (i=0; i<l; ++i)
    if (q[i] >= 63) // q30: 33 + 30
      ++(rs->q30_bases_cell_barcode);

  if (!br->n_wl) {
    barcode_tag_create (b, q, l, NULL, buf);
    return 0;
  }

  *em = 0;
  s = ss_query (br->wl, b, l, br->dist, em, thid);
  if (s == NULL)
    return -1;

  barcode_tag_create (b, q, l, s, buf);

  return 0;
}

static int
barcode_process (pe_seq_t * seq, int nbr, struct bcode_reg * brs,
    char * tag, char * raw_tag, char * raw_qual_tag, int thid,
		str_t * s, int * oem, read_stat_t * rs)
{
  int i;
  int em;
  int ret;

  str_clear (s);
  str_clear (s+1);
  str_clear (s+2);

  *oem = 1;

	ret = 0;
  if (process1barcode(seq,brs,thid,s,&em,rs) != 0)
		ret = -1;
  if (em == 0)
    *oem = 0;

  for (i=1; i<nbr; ++i) {
    if (process1barcode(seq,brs+i,thid,s,&em,rs) != 0)
			ret = -1;
    if (em == 0)
      *oem = 0;
  }

	if (ret)
		return ret;

  if (tag)
    read_name_update (seq->n, tag, s[0].s, s[0].l);

  if (raw_tag)
    read_name_update (seq->n, raw_tag, s[1].s, s[1].l);

  if (raw_qual_tag)
    read_name_update (seq->n, raw_qual_tag, s[2].s, s[2].l);

  return 0;
}

static int
seq_process (pe_seq_t * seq, struct bcode_reg * br, str_t * bases, str_t * quals, int thid, read_stat_t * rs)
{
  char * b;
  char * q;
  int i, l;
  int ridx;
  int start;
  int end;
	int nbad;

	if (br == NULL)
		return 0;

  ridx = br->rd - 1;
  start = br->start;
  end = br->end;

  b = seq->b[ridx] + start - 1;
  q = seq->q[ridx] + start - 1;
  l = end - start + 1;

  str_assign2 (bases, b, l);
  str_assign2 (quals, q, l);

  rs->bases_read += l;
  for (i=0; i<l; ++i)
    if (q[i] >= 63) // q30: 33 + 30
      ++(rs->q30_bases_read);

  if (opt->dropN) {
    for (i=0; i<l; ++i) {
      if (b[i] == 'N') {
        seq->flag = FQ_FLAG_READ_LOWQ;
        break;
      }
    }
  }

  if (seq->flag != 0)
    return 0;

	for (i=nbad=0; i<15&&i<l; ++i)
		if (q[i]<43 && (++nbad)>2) { // q10: 33 + 10
			seq->flag = FQ_FLAG_READ_LOWQ;
			return -1;
		}

	return 0;
}

static void
read_1_process (pe_seq_t * seq, struct bcode_reg * br1, struct bcode_reg * br2, int thid, read_stat_t * rs, str_t * r)
{
  int i;
  int avg;

	str_clear (r);
	str_clear (r+1);
	str_clear (r+2);
	str_clear (r+3);

  if (seq_process(seq,br1,r,r+1,thid,rs) != 0)
		return;
  if (seq_process(seq,br2,r+2,r+3,thid,rs) != 0)
		return;

  assert (r[0].l == r[1].l);
  seq->l[0] = r[0].l;
  pe_seq_resize (seq, seq->l[0], -1);
  memcpy (seq->b[0], r[0].s, seq->l[0]+1);
  memcpy (seq->q[0], r[1].s, seq->l[0]+1);

  if (r[2].l <= 0)
    return;
  assert (r[2].l == r[3].l);
  seq->l[1] = r[2].l;
  pe_seq_resize (seq, -1, seq->l[1]);
  memcpy (seq->b[1], r[2].s, seq->l[1]+1);
  memcpy (seq->q[1], r[3].s, seq->l[1]+1);

  if (seq->flag != 0)
    return;

  if (opt->min_qual>0 && seq->l[0]>0) {
    for (i=avg=0; i<seq->l[0]; ++i)
      avg += seq->q[0][i] - 33;
    if (avg/seq->l[0] < opt->min_qual) {
      seq->flag = FQ_FLAG_READ_LOWQ;
      return;
    }

    if (seq->l[1] > 0) {
      for (i=avg=0; i<seq->l[1]; ++i)
        avg += seq->q[1][i] - 33;
      if (avg/seq->l[1] < opt->min_qual) {
        seq->flag = FQ_FLAG_READ_LOWQ;
        return;
      }
    }
  }
}

static void
add1read (qual_stat_t * st, read_stat_t * rs, int thid)
{
	st->bases_cell_barcode[thid] += rs->bases_cell_barcode;
	st->q30_bases_cell_barcode[thid] += rs->q30_bases_cell_barcode;

	st->bases_umi[thid] += rs->bases_umi;
	st->q30_bases_umi[thid] += rs->q30_bases_umi;

	st->bases_read[thid] += rs->bases_read;
	st->q30_bases_read[thid] += rs->q30_bases_read;
}

static void *
read_process_core (void * d)
{
	char bc_str[256];
  int * signal;
  int em;
  int thid;
  int32_t i;
	str_t bs;
  str_t * buf;
  pe_seq_t * seq;
  xth_data_t * xd;
	read_stat_t rs;

  xd = (xth_data_t *) d;
  thid = xd->thid;
  signal = xd->signal;
	bs.s = bc_str;

	buf = (str_t *) ckmalloc (8 * sizeof(str_t));
	for (i=0; i<8; ++i)
		str_init2 (buf+i);

  for (;;) {
    if (*signal == SIG_WORK) {
			/*
			n = (blk->n+1) / opt->nth;
			start = thid * n;
			end = (thid+1) * n;
			if (end > blk->n)
				end = blk->n;
			*/
      for (i=0; i<blk->n; ++i) {
				if (i % opt->nth != thid)
					continue;
        seq = blk->seqs + i;
        seq->flag = 0;

        /*
        int zzz = 0;
        if (memcmp(seq->n->s,"@E100020122L1C001R00100401739",29) == 0)
          zzz = 1;
          */

				memset (&rs, 0, sizeof(read_stat_t));

        if (cfg->UMI)
          umi_process (seq, thid, &rs);

        if (cfg->cell_barcodes) {
          if (barcode_process (seq, cfg->n_cell_barcode, cfg->cell_barcodes,
              cfg->cell_barcode_tag, cfg->raw_cell_barcode_tag,
              cfg->raw_cell_barcode_qual_tag, thid, buf, &em, &rs) != 0) {
            seq->flag |= PE_FLAG_CBC_FAIL;
						rs.bases_cell_barcode = 0;
						rs.q30_bases_cell_barcode = 0;
						add1read (&st, &rs, thid);
            continue;
          }

          if (em == 1)
            seq->flag = PE_CBC_EXACT_MATCH;

          if (opt->cbdis_file && (!seq->flag||seq->flag==PE_CBC_EXACT_MATCH) && buf->l) {
						strcpy (bc_str, buf->s);
						bc_str[buf->l] = '\0';
						bs.l = buf->l;
          }
        }

        if (cfg->read_1)
          read_1_process (seq, cfg->read_1, cfg->read_2, thid, &rs, buf);

				if (seq->flag == FQ_FLAG_READ_LOWQ)
					continue;

				add1read (&st, &rs, thid);

				str_hash_add (bc_hashs[thid], &bs);
      }

      *signal = SIG_WAIT;
    } else if (*signal == SIG_STOP) {
      *signal = SIG_WAIT;
      break;
    }
		usleep (10);
  }

	return (void *) 0;
}

static void
check_parameters (void)
{
	assert (cfg->read_1 != NULL);

	if (cfg->cell_barcodes)
  	assert (cfg->cell_barcode_tag!=NULL || cfg->raw_cell_barcode_tag!=NULL);

  if (cfg->n_sample_barcode>0 || cfg->sample_barcodes)
    err_mesg ("Sample barcodes are not supported right now!");
}

static void
reads_output (void)
{
  int32_t i;
  pe_seq_t * seq;

  for (i=0; i<blk->n; ++i) {
    seq = blk->seqs + i;
    ++st.raw_reads;

    /*
    int zzz = 0;
    if (memcmp(seq->n->s,"@E100020122L1C001R00100401739",29) == 0)
      zzz = 1;
          */

    if (seq->flag & PE_FLAG_SBC_FAIL) {
      ++st.filter_by_sample;
      continue;
    }

    if (seq->flag & PE_FLAG_CBC_FAIL) {
      ++st.filter_by_bc;
      continue;
    }

    if (seq->flag == FQ_FLAG_READ_LOWQ) {
      ++st.filter_by_lowqual;
      continue;
    }

    ++st.reads_pass_qc;
    if (seq->flag & PE_CBC_EXACT_MATCH)
      ++st.bc_exact_match;

    if (opt->flag & OPT_FQ1_OUT_GZ)
      gzprintf (opt->gfp[0], "%s\n%s\n+\n%s\n", seq->n->s, seq->b[0], seq->q[0]);
    else
      fprintf (opt->fp[0], "%s\n%s\n+\n%s\n", seq->n->s, seq->b[0], seq->q[0]);

    if (opt->out_file[1] && seq->l[1]>0) {
      if (opt->flag & OPT_FQ2_OUT_GZ)
        gzprintf (opt->gfp[1], "%s\n%s\n+\n%s\n", seq->n->s, seq->b[1], seq->q[1]);
      else
        fprintf (opt->fp[1], "%s\n%s\n+\n%s\n", seq->n->s, seq->b[1], seq->q[1]);
    }
  }
}

static int
bc_cmp (const void * a, const void * b)
{
  xh_item_t * pa = (xh_item_t *) a;
  xh_item_t * pb = (xh_item_t *) b;

  if (pa->multi < pb->multi)
    return 1;
  else if (pa->multi > pb->multi)
    return -1;
	else
		return 0;
}

static void
cell_barcode_stat_dump (const char * file)
{
  uint64_t i;
  FILE * out;
  str_t * s;
  xh_t * h;
  xh_item_t * ptr;
  str_hash_t * hash;

  hash = bc_hashs[0];
  for (i=1; i<opt->nth; ++i)
    str_hash_merge (hash, bc_hashs[i]);

  h = hash->hash;
  qsort (h->pool, h->cnt, sizeof(xh_item_t), bc_cmp);

  out = ckopen (file, "w");
  for (i=0; i<h->cnt; ++i) {
    ptr = h->pool + i;
    s = (str_t *) ptr->key;
    fprintf (out, "%s\t%u\n", s->s, ptr->multi);
  }
  fclose (out);
}

static void
read_stat_init (void)
{
  memset (&st, 0, sizeof(qual_stat_t));

  st.bases_cell_barcode = (uint64_t *) ckalloc (opt->nth, sizeof(uint64_t));
  st.q30_bases_cell_barcode = (uint64_t *) ckalloc (opt->nth, sizeof(uint64_t));

  st.bases_sample_barcode = (uint64_t *) ckalloc (opt->nth, sizeof(uint64_t));
  st.q30_bases_sample_barcode = (uint64_t *) ckalloc (opt->nth, sizeof(uint64_t));

  st.bases_umi = (uint64_t *) ckalloc (opt->nth, sizeof(uint64_t));
  st.q30_bases_umi = (uint64_t *) ckalloc (opt->nth, sizeof(uint64_t));

  st.bases_read = (uint64_t *) ckalloc (opt->nth, sizeof(uint64_t));
  st.q30_bases_read = (uint64_t *) ckalloc (opt->nth, sizeof(uint64_t));
}

static void
report_file_create (const char * file)
{
  char out_file[4096];
  double sum;
  double q30_sum;
  FILE * out;

  if (str_cmp_tail(file,".csv",4) != 0)
    sprintf (out_file, "%s.csv", file);
  else
    strcpy (out_file, file);

  out = ckopen (out_file, "w");
  fprintf (out, "Number of Fragments,%" PRIu64 "\n", st.raw_reads);
  fprintf (out, "Fragments pass QC,%" PRIu64 "\n", st.reads_pass_qc);
  fprintf (out, "Fragments with Exactly Matched Barcodes,%" PRIu64 "\n", st.bc_exact_match);
  fprintf (out, "Fragments with Failed Barcodes,%" PRIu64 "\n", st.filter_by_bc);
  fprintf (out, "Fragments Filtered on Low Quality,%" PRIu64 "\n", st.filter_by_lowqual);
  fprintf (out, "Fragments Filtered on Unknown Sample Barcodes,%" PRIu64 "\n", st.filter_by_sample);

  q30_sum = (double) sum (u64, st.q30_bases_cell_barcode, opt->nth);
  sum = (double) sum (u64, st.bases_cell_barcode, opt->nth);
  fprintf (out, "Q30 bases in Cell Barcode,%.1f%%\n", sum<1e-2 ? 0 : 100*q30_sum/sum);

  q30_sum = (double) sum (u64, st.q30_bases_sample_barcode, opt->nth);
  sum = (double) sum (u64, st.bases_sample_barcode, opt->nth);
  fprintf (out, "Q30 bases in Sample Barcode,%.1f%%\n", sum<1e-2 ? 0 : 100*q30_sum/sum);

  q30_sum = (double) sum (u64, st.q30_bases_umi, opt->nth);
  sum = (double) sum (u64, st.bases_umi, opt->nth);
  fprintf (out, "Q30 bases in UMI,%.1f%%\n", sum<1e-2 ? 0 : 100*q30_sum/sum);

  q30_sum = (double) sum (u64, st.q30_bases_read, opt->nth);
  sum = (double) sum (u64, st.bases_read, opt->nth);
  fprintf (out, "Q30 bases in Reads,%.1f%%\n", sum<1e-2 ? 0 : 100*q30_sum/sum);
  fclose (out);
}

int
fastq_prase_barcodes (int argc, char * argv[])
{
  int i;
	int64_t n_reads;
  time_t time_start;
  xth_engine_t * engine;

  time (&time_start);

  opt = parse_options (argc, argv);
  cfg = config_init2 (opt->cfg_file, opt->nth);

  check_parameters ();

  assert (opt->nth > 0);

	/*
  cbc_buf = (str_t *) ckmalloc (3 * opt->nth * sizeof(str_t));
  read_buf = (str_t *) ckmalloc (4 * opt->nth * sizeof(str_t));
	for (i=0; i<3*opt->nth; ++i)
		str_init2 (cbc_buf+i);
  for (i=0; i<4*opt->nth; ++i)
    str_init2 (read_buf+i);
		*/

  if (opt->cbdis_file) {
    bc_hashs = (str_hash_t **) ckmalloc (opt->nth * sizeof(str_hash_t *));
    for (i=0; i<opt->nth; ++i)
      bc_hashs[i] = str_hash_init ();
  }

  read_stat_init ();

  engine = xth_engine_start (read_process_core, NULL, opt->nth, SIG_WAIT);

	n_reads = 0;
  blk = pefq_blk_init (BATCH_RD_CNT);
  for (i=0; i<opt->n; ++i) {
    reader = pefq_reader_init (opt->in_file[0][i], opt->in_file[1][i]);
    while (pefq_reader_load_block (reader, blk) == 0) {
      cnter = 0;

			//time (&beg);
      xth_engine_send_signal (engine, SIG_WORK, SIG_WAIT);
			//printf ("Fastq processing cost: %lds\n", time(NULL)-beg);

			//time (&beg);
      reads_output ();
			//printf ("Fastq dumping cost: %lds\n", time(NULL)-beg);

			n_reads += blk->n;
			printf ("%ld reads have been processed\n", n_reads);
    }
  }

  xth_engine_stop (engine, SIG_STOP, SIG_WAIT);

  if (opt->cbdis_file)
    cell_barcode_stat_dump (opt->cbdis_file);

  if (opt->report_file)
    report_file_create (opt->report_file);

  if (opt->flag & OPT_FQ1_OUT_GZ)
    gzclose (opt->gfp[0]);
  else
    fclose (opt->fp[0]);

  if (opt->out_file[1]) {
    if (opt->flag & OPT_FQ2_OUT_GZ)
      gzclose (opt->gfp[1]);
    else
      fclose (opt->fp[1]);
  }

  config_destroy2 (cfg);

  printf ("Program cost: %lds\n", time(NULL)-time_start);

  return 0;
}
