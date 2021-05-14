/***********************************************************
  *File Name: 
  *Description: 
  *Author: Chen Xi
  *Email: chenxi1@genomics.cn
  *Create Time: 2019-09-02 16:00:51
  *Edit History: 
***********************************************************/

#ifndef XDK_MBASIC_H
#define XDK_MBASIC_H

#include <stdint.h>

#define NEGA_INF (-INFINITY)

#define MATH_BASIC(name, type_t) \
  static inline type_t sum_##name (type_t * arr, int32_t n) { \
    int32_t i; \
    type_t sum = 0; \
    for (i=0; i<n; ++i) \
      sum += arr[i]; \
    return sum; \
  } \
  \
  static inline void zero_##name (type_t * arr, int32_t n) { \
    memset (arr, 0, n*sizeof(type_t)); \
  }

#define sum(name,arr,n) sum_##name((arr),(n))
#define zero(name,arr,n) zero_##name((arr),(n))

MATH_BASIC (i32, int32_t);
MATH_BASIC (u32, uint32_t);
MATH_BASIC (i64, int64_t);
MATH_BASIC (u64, uint64_t);
MATH_BASIC (dbl, double);

static inline int32_t
fast_round (double d)
{
  return (d>0.0) ? (int32_t)(d+0.5) : (int32_t)(d-0.5);
}

#ifdef __cplusplus
extern "C" {
#endif

  double digamma (double x);
  double log_gamma (double x);

  // add two arrays element-by-element
  void ebe_add_d (double * sum, double * a, double * b, int32_t n);

  double log_sum_exp (double * log_values, int32_t n);

  double distance1_d (double * v1, double * v2, int32_t n);

  double log_dirichlet_normalization (double * a, int32_t n);

  double log1p (double d);

  double ln_to_log10 (double d);
  double log10_to_ln (double d);

  double fast_bernoulli_entropy (double d);

#ifdef __cplusplus
}
#endif

/**
 * Jacobian log identiti for differences up to MAX_TOLERANCE
 */

struct jac_log_tbl_s;
typedef struct jac_log_tbl_s jac_log_tbl_t;

struct jac_log_tbl_s {
  double max_tolerance;
  double table_step;
  double inv_step;
  double * cache;
  int32_t l_cache;
};

#ifdef __cplusplus
extern "C" {
#endif

  jac_log_tbl_t * jac_log_tbl_init (void);

  void jac_log_tbl_free (jac_log_tbl_t * t);

  double jac_log_tbl_get (jac_log_tbl_t * t, double diff);

#ifdef __cplusplus
}
#endif

/**
 * gaussian kernel
 */

struct gaussian_kernel_s;
typedef struct gaussian_kernel_s gaussian_kernel_t;

#ifdef __cplusplus
extern "C" {
#endif

  gaussian_kernel_t * gaussian_kernel_init (int32_t mean, double sigma);

  double gaussian_kernal_get (gaussian_kernel_t * gk, int32_t pos);

#ifdef __cplusplus
}
#endif

#endif
