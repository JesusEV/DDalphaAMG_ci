/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#ifndef MAIN_POST_DEF_PRECISION_HEADER
  #define MAIN_POST_DEF_PRECISION_HEADER
  
  #include "coarse_oddeven_PRECISION.h"
  #include "dirac_PRECISION.h"
  #include "coarse_operator_PRECISION.h"
  #include "block_jacobi_PRECISION.h"

  static inline void apply_operator_PRECISION( vector_PRECISION output, vector_PRECISION input, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

#ifdef BLOCK_JACOBI
    if ( l->level==0 && p->block_jacobi_PRECISION.BJ_usable==1 ) {
      p->eval_operator( l->p_PRECISION.block_jacobi_PRECISION.xtmp, input, p->op, l, threading );
      block_jacobi_apply_PRECISION( output, l->p_PRECISION.block_jacobi_PRECISION.xtmp, p, l, threading );

      //printf("BJ usable = %d\n", p->block_jacobi_PRECISION.BJ_usable);

      //START_MASTER(threading)
      //printf0("SHAIT1\n");
      //printf0("BJ usable = %d\n", p->block_jacobi_PRECISION.BJ_usable);
      //END_MASTER(threading)

      // --------------------------------------------------------------------------------
      //START_MASTER(threading)
      /*
      if ( p->block_jacobi_PRECISION.BJ_usable==1 ) {
        {
          PRECISION tmpx1, tmpx2;

          //int start = p->v_start;
          //int end = p->v_end;

          int start, end;
          compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

          //local_gmres_PRECISION_struct *loc_p = &(p->block_jacobi_PRECISION.local_p);

          int size = p->v_end - p->v_start;
          int exp_fctr = 10;

          START_MASTER(threading)
          p->block_jacobi_PRECISION.xxxtmp[0] = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );
          p->block_jacobi_PRECISION.xxxtmp[1] = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );
          p->block_jacobi_PRECISION.xxxtmp[2] = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );
          END_MASTER(threading)

          SYNC_MASTER_TO_ALL(threading)
          SYNC_CORES(threading)

          vector_PRECISION solution = p->block_jacobi_PRECISION.xxxtmp[0];
          vector_PRECISION rhs = p->block_jacobi_PRECISION.xxxtmp[1];
          vector_PRECISION x = p->block_jacobi_PRECISION.xxxtmp[2];

          SYNC_MASTER_TO_ALL(threading)
          SYNC_CORES(threading)

          START_MASTER(threading)
          vector_PRECISION_define_random( solution, p->v_start, p->v_end, l );
          END_MASTER(threading)

          START_MASTER(threading)
          printf0("solution[0] = %f+i%f\n", creal(solution[0]), cimag(solution[0]));
          END_MASTER(threading)

          SYNC_MASTER_TO_ALL(threading)
          SYNC_CORES(threading)

          START_MASTER(threading)
          printf0("SHAIT2\n");
          END_MASTER(threading)

          p->eval_operator(rhs, solution, p->op, l, threading);

          START_MASTER(threading)
          printf0("SHAIT2\n");
          END_MASTER(threading)

          SYNC_MASTER_TO_ALL(threading)
          SYNC_CORES(threading)

          // x ~ w
          //local_apply_polyprec_PRECISION( x, NULL, rhs, 0, l, threading );
          block_jacobi_apply_PRECISION( x, rhs, p, l, threading );

          START_MASTER(threading)
          printf0("rhs[0] = %f+i%f\n", creal(rhs[0]), cimag(rhs[0]));
          END_MASTER(threading)

          START_MASTER(threading)
          printf0("x[0] = %f+i%f\n", creal(x[0]), cimag(x[0]));
          END_MASTER(threading)

          START_MASTER(threading)
          p->block_jacobi_PRECISION.xxxtmp[3] = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );
          END_MASTER(threading)

          SYNC_MASTER_TO_ALL(threading)
          SYNC_CORES(threading)

          vector_PRECISION diff_sol = p->block_jacobi_PRECISION.xxxtmp[3];

          SYNC_MASTER_TO_ALL(threading)
          SYNC_CORES(threading)

          //vector_PRECISION diff_sol = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );

          vector_PRECISION_minus( diff_sol, x, solution, start, end, l );

          START_MASTER(threading)
          printf0("diff_sol[0] = %f+i%f\n", creal(diff_sol[0]), cimag(diff_sol[0]));
          END_MASTER(threading)

          tmpx1 = global_norm_PRECISION( diff_sol, start, end, l, threading );
          tmpx2 = global_norm_PRECISION( solution, start, end, l, threading );

          printf0("g (proc=%d) ---> approx rel error BJ = %f\n", g.my_rank, tmpx1/tmpx2);

          START_MASTER(threading)
          free(solution);
          free(rhs);
          free(x);
          free(diff_sol);
          END_MASTER(threading)
        }
      }
      //END_MASTER(threading)

      START_MASTER(threading)
      printf0("SHAIT3\n");
      END_MASTER(threading)

      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)

      START_MASTER(threading)
      MPI_Barrier( MPI_COMM_WORLD );
      MPI_Finalize();
      END_MASTER(threading)

      SYNC_MASTER_TO_ALL(threading)
      SYNC_CORES(threading)

      exit(0);
      */
      // --------------------------------------------------------------------------------
    } else {
      p->eval_operator( output, input, p->op, l, threading );
    }
#else
    p->eval_operator( output, input, p->op, l, threading );
#endif
  }
  
  static inline void apply_operator_dagger_PRECISION( vector_PRECISION output, vector_PRECISION input, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      tau1_gamma5_PRECISION( l->vbuf_PRECISION[6], input, l, threading );
    } else
#endif
      {
        gamma5_PRECISION( l->vbuf_PRECISION[6], input, l, threading );
#ifdef HAVE_TM
        //TODO: change_mu_sign_PRECISION( p->op, l, threading );
#endif
      }

    apply_operator_PRECISION( l->vbuf_PRECISION[7], l->vbuf_PRECISION[6], p, l, threading );

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      tau1_gamma5_PRECISION( output, l->vbuf_PRECISION[7], l, threading );
    } else
#endif
      {
        gamma5_PRECISION( output, l->vbuf_PRECISION[7], l, threading );
#ifdef HAVE_TM
        //TODO: change_mu_sign_PRECISION( p->op, l, threading );
#endif
      }
    
  }

  static inline void test0_PRECISION( char* format, int depth, PRECISION test ) {
    if ( g.my_rank == 0 && g.print >= 0 ) {
      if ( test > EPS_PRECISION )
        printf("\x1b[31m");
      printf(format, depth, test);
      if ( test > EPS_PRECISION )
        printf("\x1b[0m");
      if ( test > g.test )
        g.test = test;
      fflush(0);
    }
  }
  
#endif
