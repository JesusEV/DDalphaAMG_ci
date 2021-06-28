// TODO : include in main.h : sparse_mv_float.h and sparse_mv_double.h

#ifndef SPARSE_MV_PRECISION_HEADER
  #define SPARSE_MV_PRECISION_HEADER

  struct Thread;

  void spmv_PRECISION( vector_PRECISION out, vector_PRECISION in, vector_PRECISION A, int *Is, int *Js,
                       int n, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading );

#endif
