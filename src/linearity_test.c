#include <float.h>
#include "main.h"

// own random generator (xorshift64*), so the global rand() sequence and thus
// any random right-hand side of the following solve stays unchanged
static double lt_rand( unsigned long long *s ) {
  *s ^= *s >> 12; *s ^= *s << 25; *s ^= *s >> 27;
  return (double)( (*s * 2685821657736338717ULL) >> 11 ) * ( 1.0/9007199254740992.0 ) - 0.5;
}

static void lt_fill_random( vector_double v, long int n, unsigned long long *s ) {
  for ( long int i=0; i<n; i++ ) {
    double re = lt_rand( s );
    double im = lt_rand( s );
    v[i] = re + im*_Complex_I;
  }
}

// applies exactly the preconditioner the outer solver uses, with zero initial guess
static void lt_apply_prec( vector_double phi, vector_double eta, level_struct *l, struct Thread *threading ) {
  if ( g.mixed_precision == 2 ) {
    // fgmres_MP: preconditioner is vcycle_float on the translated float vectors
    trans_float( l->sbuf_float[0], eta, l->s_float.op.translation_table, l, threading );
    vcycle_float( l->sbuf_float[1], NULL, l->sbuf_float[0], _NO_RES, l, threading );
    trans_back_float( phi, l->sbuf_float[1], l->s_float.op.translation_table, l, threading );
  } else {
    // fgmres_double: preconditioner() (handles mixed_precision 0 and 1)
    preconditioner( phi, NULL, eta, _NO_RES, l, threading );
  }
  SYNC_MASTER_TO_ALL(threading)
}

// returns || a - b || / || b ||, d is workspace
static double lt_rel_diff( vector_double a, vector_double b, vector_double d,
                           int start, int end, level_struct *l, struct Thread *threading ) {
  vector_double_minus( d, a, b, start, end, l );
  SYNC_MASTER_TO_ALL(threading)
  double nd = global_norm_double( d, 0, l->inner_vector_size, l, threading );
  double nb = global_norm_double( b, 0, l->inner_vector_size, l, threading );
  return ( nb > 0 ) ? nd/nb : nd;
}


void linearity_test( level_struct *l, struct Thread *threading ) {

  if ( g.method <= 0 ) {
    START_MASTER(threading)
    printf0("linearity test: skipped, method %d uses no multigrid preconditioner\n", g.method );
    END_MASTER(threading)
    return;
  }

  int start = threading->start_index[l->depth];
  int end   = threading->end_index[l->depth];
  long int vs = l->inner_vector_size;
#ifdef HAVE_TM1p1
  vs *= 2;
#endif
  const complex_double alpha = 0.314 + 1.592*_Complex_I;

  // one contiguous block for all 11 vectors
  vector_double buf = NULL;
  PUBLIC_MALLOC( buf, complex_double, 11*vs );
  vector_double x    = buf,      y    = buf+  vs, ax   = buf+2*vs, xpy = buf+3*vs,
                Mx   = buf+4*vs, My   = buf+5*vs, Max  = buf+6*vs, Mxpy = buf+7*vs,
                Mx2  = buf+8*vs, ref  = buf+9*vs, d    = buf+10*vs;

  // save timers / counters that the preconditioner calls would otherwise falsify
  double save_coarse_time = 0, save_coarsest_time = 0;
#if defined(MUMPS_ADDS) || defined(COARSE_SCALAP)
  double save_cs_time = 0; int save_cs_number = 0;
#endif
#ifdef COARSE_SCALAP
  double save_comm_time = 0;
#endif

  START_LOCKED_MASTER(threading)
  save_coarse_time = g.coarse_time; save_coarsest_time = g.coarsest_time;
#if defined(MUMPS_ADDS) || defined(COARSE_SCALAP)
  save_cs_time = g.coarsest_solve_time; save_cs_number = g.coarsest_solve_number;
#endif
#ifdef COARSE_SCALAP
  save_comm_time = g.scalap_comm_time;
#endif
  // different random vectors on every rank, but reproducible from run to run
  unsigned long long seed = 0x9E3779B97F4A7C15ULL ^ ( 1000003ULL * (unsigned long long)(g.my_rank+1) );
  lt_fill_random( x, vs, &seed );
  lt_fill_random( y, vs, &seed );
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // test inputs
  vector_double_scale( ax, x, alpha, start, end, l );  // ax  = alpha*x
  vector_double_plus( xpy, x, y, start, end, l );      // xpy = x + y
  SYNC_MASTER_TO_ALL(threading)

  double t0 = MPI_Wtime();
  lt_apply_prec( Mx,   x,   l, threading );   // first call
  lt_apply_prec( My,   y,   l, threading );
  lt_apply_prec( Max,  ax,  l, threading );
  lt_apply_prec( Mxpy, xpy, l, threading );
  lt_apply_prec( Mx2,  x,   l, threading );   // last call, same input as the first one
  double t1 = MPI_Wtime();

  double nx  = global_norm_double( x,  0, l->inner_vector_size, l, threading );
  double nMx = global_norm_double( Mx, 0, l->inner_vector_size, l, threading );

  vector_double_scale( ref, Mx, alpha, start, end, l ); // ref = alpha*M(x)
  SYNC_MASTER_TO_ALL(threading)
  double e_hom = lt_rel_diff( Max, ref, d, start, end, l, threading );

  vector_double_plus( ref, Mx, My, start, end, l );     // ref = M(x) + M(y)
  SYNC_MASTER_TO_ALL(threading)
  double e_add = lt_rel_diff( Mxpy, ref, d, start, end, l, threading );

  double e_rep = lt_rel_diff( Mx2, Mx, d, start, end, l, threading );

  // heuristic threshold: a few orders of magnitude above the unit roundoff of
  // the precision the preconditioner runs in
  double eps = g.mixed_precision ? FLT_EPSILON : DBL_EPSILON;
  double thr = 1E3*eps;
  double e_max = MAX( e_hom, MAX( e_add, e_rep ) );

  START_LOCKED_MASTER(threading)
  g.coarse_time = save_coarse_time; g.coarsest_time = save_coarsest_time;
#if defined(MUMPS_ADDS) || defined(COARSE_SCALAP)
  g.coarsest_solve_time = save_cs_time; g.coarsest_solve_number = save_cs_number;
#endif
#ifdef COARSE_SCALAP
  g.scalap_comm_time = save_comm_time;
#endif

  printf0("\n+----------------------------------------------------------+\n");
  printf0("|           linearity test of the preconditioner           |\n");
  printf0("+----------------------------------------------------------+\n");
  printf0("| alpha = %+.2f%+.2fi, 5 applications, %8.3f s           |\n", creal(alpha), cimag(alpha), t1-t0 );
  printf0("| precision of M: %-6s  unit roundoff: %e      |\n", g.mixed_precision?"single":"double", eps );
  printf0("| ||M(x)|| / ||x||                         : %e |\n", nMx/nx );
  printf0("| homogeneity   ||M(ax)-aM(x)||/||aM(x)||  : %e |\n", e_hom );
  printf0("| additivity    ||M(x+y)-Mx-My||/||Mx+My|| : %e |\n", e_add );
  printf0("| repeatability ||M(x)_5-M(x)_1||/||M(x)|| : %e |\n", e_rep );
  printf0("+----------------------------------------------------------+\n");
  if ( e_max < thr )
    printf0("| => M is linear and stationary up to rounding (< %.1e) |\n", thr );
  else
    printf0("| => M is NOT linear/stationary (max %.2e > %.1e)  |\n", e_max, thr );
  printf0("+----------------------------------------------------------+\n\n");
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  PUBLIC_FREE( buf, complex_double, 11*vs );
}
