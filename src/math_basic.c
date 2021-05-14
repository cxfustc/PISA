/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-09-02 16:01:23
  *Edit History: 
***********************************************************/

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "util.h"
#include "math_basic.h"

static int LANCZOS_LENGTH = 15;
static double LANCZOS[15] = {
  0.99999999999999709182,
  57.156235665862923517,
  -59.597960355475491248,
  14.136097974741747174,
  -0.49191381609762019978,
  .33994649984811888699e-4,
  .46523628927048575665e-4,
  -.98374475304879564677e-4,
  .15808870322491248884e-3,
  -.21026444172410488319e-3,
  .21743961811521264320e-3,
  -.16431810653676389022e-3,
  .84418223983852743293e-4,
  -.26190838401581408670e-4,
  .36899182659531622704e-5,
};

static double LANCZOS_G = 607.0 / 128.0;
static double HALF_LOG_2_PI = 0.9189385; // 0.5 * log(2*pi)

static double S_LIMIT = 1.0e-5;
static double C_LIMIT = 49.0;
static double GAMMA   = 0.577215664901532860606512090082;

/*
 * Constants for the computation of double invGamma1pm1(double).
 * Copied from DGAM1 in the NSWC library.
 */

/** The constant {@code A0} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_A0 = .611609510448141581788E-08;

/** The constant {@code A1} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_A1 = .624730830116465516210E-08;

/** The constant {@code B1} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_B1 = .203610414066806987300E+00;

/** The constant {@code B2} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_B2 = .266205348428949217746E-01;

/** The constant {@code B3} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_B3 = .493944979382446875238E-03;

/** The constant {@code B4} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_B4 = -.851419432440314906588E-05;

/** The constant {@code B5} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_B5 = -.643045481779353022248E-05;

/** The constant {@code B6} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_B6 = .992641840672773722196E-06;

/** The constant {@code B7} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_B7 = -.607761895722825260739E-07;

/** The constant {@code B8} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_B8 = .195755836614639731882E-09;

/** The constant {@code P0} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_P0 = .6116095104481415817861E-08;

/** The constant {@code P1} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_P1 = .6871674113067198736152E-08;

/** The constant {@code P2} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_P2 = .6820161668496170657918E-09;

/** The constant {@code P3} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_P3 = .4686843322948848031080E-10;

/** The constant {@code P4} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_P4 = .1572833027710446286995E-11;

/** The constant {@code P5} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_P5 = -.1249441572276366213222E-12;

/** The constant {@code P6} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_P6 = .4343529937408594255178E-14;

/** The constant {@code Q1} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_Q1 = .3056961078365221025009E+00;

/** The constant {@code Q2} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_Q2 = .5464213086042296536016E-01;

/** The constant {@code Q3} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_Q3 = .4956830093825887312020E-02;

/** The constant {@code Q4} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_Q4 = .2692369466186361192876E-03;

/** The constant {@code C} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C = -.422784335098467139393487909917598E+00;

/** The constant {@code C0} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C0 = .577215664901532860606512090082402E+00;

/** The constant {@code C1} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C1 = -.655878071520253881077019515145390E+00;

/** The constant {@code C2} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C2 = -.420026350340952355290039348754298E-01;

/** The constant {@code C3} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C3 = .166538611382291489501700795102105E+00;

/** The constant {@code C4} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C4 = -.421977345555443367482083012891874E-01;

/** The constant {@code C5} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C5 = -.962197152787697356211492167234820E-02;

/** The constant {@code C6} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C6 = .721894324666309954239501034044657E-02;

/** The constant {@code C7} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C7 = -.116516759185906511211397108401839E-02;

/** The constant {@code C8} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C8 = -.215241674114950972815729963053648E-03;

/** The constant {@code C9} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C9 = .128050282388116186153198626328164E-03;

/** The constant {@code C10} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C10 = -.201348547807882386556893914210218E-04;

/** The constant {@code C11} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C11 = -.125049348214267065734535947383309E-05;

/** The constant {@code C12} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C12 = .113302723198169588237412962033074E-05;

/** The constant {@code C13} defined in {@code DGAM1}. */
static double INV_GAMMA1P_M1_C13 = -.205633841697760710345015413002057E-06;

static inline double
digamma4 (double x)
{
  // use method 4 (acurate to O(1/x^8))
  double inv = 1 / (x * x);
  //            1      1         1         1
  // log(x) -  --- - ------ + ------- - -------
  //           2 x   12 x^2   120 x^4   252 x^6
  return log(x) - 0.5 / x - inv * ((1.0 / 12) - inv * (1.0 / 120 - inv / 252));
}

/**
 * Returns the Lanczos approximation used to compute the gamma function.
 * The Lancos approximation is related to the Gamma function by the
 * following equation.
 *
 * gamma(x) = sqrt(x*pi) / x * (x+g+0.5) ^ (x+0.5)
 *            * exp(-x-g-0.5) * lanczos(x)
 * where g is the Lanzcos constant
 */

static double
lanczos (double x)
{
  int i;
  double sum = 0;

  for (i=LANCZOS_LENGTH-1; i>0; --i)
    sum += LANCZOS[i] / (x + i);

  return sum;
}

/**
 * Computes the value of 1 / Gamma(1+x)-1 for -0.5<=x<=1.5
 */

static double
invGamma1pm1 (double x)
{
  double ret;
  double t = x <= 0.5 ? x : (x - 0.5) - 0.5;
  double a,b,c,p,q;

  if (t < 0.0) {
    a = INV_GAMMA1P_M1_A0 + t * INV_GAMMA1P_M1_A1;
    b = INV_GAMMA1P_M1_B8;
    b = INV_GAMMA1P_M1_B7 + t * b;
    b = INV_GAMMA1P_M1_B6 + t * b;
    b = INV_GAMMA1P_M1_B5 + t * b;
    b = INV_GAMMA1P_M1_B4 + t * b;
    b = INV_GAMMA1P_M1_B3 + t * b;
    b = INV_GAMMA1P_M1_B2 + t * b;
    b = INV_GAMMA1P_M1_B1 + t * b;
    b = 1.0 + t * b;

    c = INV_GAMMA1P_M1_C13 + t * (a / b);
    c = INV_GAMMA1P_M1_C12 + t * c;
    c = INV_GAMMA1P_M1_C11 + t * c;
    c = INV_GAMMA1P_M1_C10 + t * c;
    c = INV_GAMMA1P_M1_C9 + t * c;
    c = INV_GAMMA1P_M1_C8 + t * c;
    c = INV_GAMMA1P_M1_C7 + t * c;
    c = INV_GAMMA1P_M1_C6 + t * c;
    c = INV_GAMMA1P_M1_C5 + t * c;
    c = INV_GAMMA1P_M1_C4 + t * c;
    c = INV_GAMMA1P_M1_C3 + t * c;
    c = INV_GAMMA1P_M1_C2 + t * c;
    c = INV_GAMMA1P_M1_C1 + t * c;
    c = INV_GAMMA1P_M1_C + t * c;
    if (x > 0.5) {
      ret = t * c / x;
    } else {
      ret = x * ((c + 0.5) + 0.5);
    }
  } else {
    p = INV_GAMMA1P_M1_P6;
    p = INV_GAMMA1P_M1_P5 + t * p;
    p = INV_GAMMA1P_M1_P4 + t * p;
    p = INV_GAMMA1P_M1_P3 + t * p;
    p = INV_GAMMA1P_M1_P2 + t * p;
    p = INV_GAMMA1P_M1_P1 + t * p;
    p = INV_GAMMA1P_M1_P0 + t * p;

    q = INV_GAMMA1P_M1_Q4;
    q = INV_GAMMA1P_M1_Q3 + t * q;
    q = INV_GAMMA1P_M1_Q2 + t * q;
    q = INV_GAMMA1P_M1_Q1 + t * q;
    q = 1.0 + t * q;

    c = INV_GAMMA1P_M1_C13 + (p / q) * t;
    c = INV_GAMMA1P_M1_C12 + t * c;
    c = INV_GAMMA1P_M1_C11 + t * c;
    c = INV_GAMMA1P_M1_C10 + t * c;
    c = INV_GAMMA1P_M1_C9 + t * c;
    c = INV_GAMMA1P_M1_C8 + t * c;
    c = INV_GAMMA1P_M1_C7 + t * c;
    c = INV_GAMMA1P_M1_C6 + t * c;
    c = INV_GAMMA1P_M1_C5 + t * c;
    c = INV_GAMMA1P_M1_C4 + t * c;
    c = INV_GAMMA1P_M1_C3 + t * c;
    c = INV_GAMMA1P_M1_C2 + t * c;
    c = INV_GAMMA1P_M1_C1 + t * c;
    c = INV_GAMMA1P_M1_C0 + t * c;

    if (x > 0.5) {
      ret = (t / x) * ((c - 0.5) - 0.5);
    } else {
      ret = x * c;
    }
  }

  return ret;
}

/**
 * Returns the value of log Gamma (1+x) for -0.5<=x<=1.5
 */

static double
log_gamma1p (double x)
{
  if (x < -0.5 || x > 1.5)
    err_mesg ("[log_gamma1p] -0.5<=x<=1.5!");

  return -log(invGamma1pm1(x) + 1);
}

/**
 * Computes the digamma function of x
 */

double
digamma (double x)
{
  if (x>0 && x<=S_LIMIT) {
    // use method 5from Bernardo AS103
    // acurate to O(x)
    return -GAMMA - 1/x;
  }

  if (x > C_LIMIT) {
    // use method 4 (acurate to O(1/x^8))
    return digamma4(x);
  }

  double y = x + floor(C_LIMIT-x);
  double v = digamma4 (y);
  y-=1.0; x-=1e-3;
  while (y > x) {
    v -= 1.0 / y;
    y -= 1.0;
  }

  return v;
}

/**
 * Computes the value of log Gamma of x
 */

double
log_gamma (double x)
{
  int i, n;
  double prod;
  double sum;
  double tmp;

  if (x <= 0)
    err_mesg ("[log_gamma] x must > 0!");

  if (x < 0.5)
    return log_gamma1p(x) - log(x);

  if (x <= 2.5)
    return log_gamma1p((x-0.5)-0.5);

  if (x <= 8.0) {
    n = (int) floor (x-1.5);
    prod = 1.0;
    for (i=1; i<=n; ++i)
      prod *= (x - i);
    return log_gamma1p(x-(n+1)) + log(prod);
  }

  sum = lanczos (x);
  tmp = x + LANCZOS_G + .5;

  return ((x+.5) * log(tmp)) - tmp + HALF_LOG_2_PI + log(sum/x);
}

/**
 * Adds two arrays element-by-element
 */

void
ebe_add_d (double * sum, double * a, double * b, int n_elems)
{
  for (--n_elems; n_elems>=0; --n_elems)
    sum[n_elems] = a[n_elems] + b[n_elems];
}

/**
 * Computes log(sum_i e^{a_i}) trying to avoid underflow issues by using the log-sum-exp trick
 */

double
log_sum_exp (double * log_values, int32_t n)
{
  int32_t i;
  int32_t max_idx;
  double sum;
  double max_val;

  max_idx = 0;
  max_val = log_values[0];
  for (i=1; i<n; ++i) {
    if (log_values[i] > max_val) {
      max_idx = i;
      max_val = log_values[i];
    }
  }
  if (max_val == -DBL_MAX)
    return -DBL_MAX;

  sum = 1.0;
  for (i=0; i<n; ++i) {
    if (i == max_idx)
      continue;
    sum += exp (log_values[i]-max_val);
  }

  return max_val + log(sum);
}

/**
 * Calculates the L1 (sum of abs) distance between two points
 */

double
distance1_d (double * v1, double * v2, int32_t n)
{
  int32_t i;
  double sum;

  sum = 0;
  for (i=0; i<n; ++i)
    sum += fabs (v1[i]-v2[i]);

  return sum;
}

/**
 *
 */

double
log_dirichlet_normalization (double * a, int32_t n)
{
  int32_t i;
  double log_numerator;
  double log_denominator;

  log_numerator = log_gamma (sum(dbl,a,n));
  log_denominator = 0;
  for (i=0; i<n; ++i)
    log_denominator += log_gamma (a[i]);

  return log_numerator - log_denominator;
}

double
log1p (double d)
{
  if (d < -1)
    return NAN;
  if (d == -1)
    return NEGA_INF;
  if (d == INFINITY)
    return INFINITY;

  return log(1+d);
}

static double LN10 = 2.302585;

double
ln_to_log10 (double d)
{
  return d / LN10;
}

double
log10_to_ln (double d)
{
  return d * LN10;
}

/**
 * Computes the entropy -p*ln(p) - (1-p)*ln(1-p) of a Bernoulli distribution with success probability p
 * using an extremely fast Pade approximation that is very accurate for all values of 0 <= p <= 1
 * See http://www.nezumi.demon.co.uk/consult/logx.htm
 */

double
fast_bernoulli_entropy (double d)
{
  double prod = d * (1-d);
  return prod * (11+33*prod) / (2+20*prod);
}

/**
 * Jacobian log identiti for differences up to MAX_TOLERANCE
 */

jac_log_tbl_t *
jac_log_tbl_init (void)
{
  int32_t i;
  jac_log_tbl_t * t;

  t = (jac_log_tbl_t *) ckalloc (1, sizeof(jac_log_tbl_t));
  t->max_tolerance = 8.0;
  t->table_step = 0.0001;
  t->inv_step = 1.0 / t->table_step;
  t->l_cache = (int32_t)(t->max_tolerance/t->table_step) + 1;
  t->cache = (double *) ckmalloc (t->l_cache * sizeof(double));
  for (i=0; i<t->l_cache; ++i)
    t->cache[i] = log10 (1.0 + pow(10.0,-i*t->table_step));

  return t;
}

void
jac_log_tbl_free (jac_log_tbl_t * t)
{
  free (t->cache);
  free (t);
}

double
jac_log_tbl_get (jac_log_tbl_t * t, double diff)
{
  int32_t index = fast_round (diff * t->inv_step);
  return t->cache[index];
}

/**
 * gaussian kernel
 */

struct gaussian_kernel_s {
  int32_t mean;
  int32_t band_size;
  double sigma;
  double * vals;
};

static double INV_SQRT_2_PI = 0.3989423;

static double
normal_distribution (double mean, double sd, double x)
{
  return (INV_SQRT_2_PI / sd) * exp (-(x-mean)*(x-mean) / (2.0*sd*sd));
}

gaussian_kernel_t *
gaussian_kernel_init (int32_t mean, double sigma)
{
  int32_t i;
  double sum;
  gaussian_kernel_t * gk;

  gk = (gaussian_kernel_t *) ckmalloc (sizeof(gaussian_kernel_t));
  gk->mean = mean;
  gk->sigma = sigma;
  gk->band_size = 2*mean + 1;
  gk->vals = (double *) ckmalloc (gk->band_size*sizeof(double));
  for (i=0; i<gk->band_size; ++i)
    gk->vals[i] = normal_distribution (mean, sigma, i);

  sum = 0;
  for (i=0; i<gk->band_size; ++i)
    sum += gk->vals[i];
  assert (sum != 0);
  for (i=0; i<gk->band_size; ++i)
    gk->vals[i] = gk->vals[i] / sum;

  return gk;
}

double
gaussian_kernal_get (gaussian_kernel_t * gk, int32_t pos)
{
  assert (pos>=0 && pos<gk->band_size);
  return gk->vals[pos];
}
