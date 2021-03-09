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
      //p->eval_operator( output, l->p_PRECISION.block_jacobi_PRECISION.xtmp, p->op, l, threading );
      //block_jacobi_apply_PRECISION( l->p_PRECISION.block_jacobi_PRECISION.xtmp, input, p, l, threading );

      p->eval_operator( l->p_PRECISION.block_jacobi_PRECISION.xtmp, input, p->op, l, threading );
      block_jacobi_apply_PRECISION( output, l->p_PRECISION.block_jacobi_PRECISION.xtmp, p, l, threading );

      // --------------------------------------------------------------------------------
      /*
      if ( p->block_jacobi_PRECISION.BJ_usable==1 ) {
        {
          PRECISION tmpx1, tmpx2;
          int start = p->v_start;
          int end = p->v_end;
          //local_gmres_PRECISION_struct *loc_p = &(p->block_jacobi_PRECISION.local_p);

          int size = end-start;
          int exp_fctr = 10;

          vector_PRECISION solution = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );
          vector_PRECISION rhs = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );
          vector_PRECISION x = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );

          vector_PRECISION_define_random( solution, start, end, l );

          p->eval_operator(rhs, solution, p->op, l, threading);

          // x ~ w
          //local_apply_polyprec_PRECISION( x, NULL, rhs, 0, l, threading );
          block_jacobi_apply_PRECISION( x, rhs, p, l, threading );

          vector_PRECISION diff_sol = (vector_PRECISION) malloc( exp_fctr*size*size*sizeof(complex_PRECISION) );

          vector_PRECISION_minus( diff_sol, x, solution, start, end, l );

          tmpx1 = global_norm_PRECISION( diff_sol, p->v_start, p->v_end, l, threading );
          tmpx2 = global_norm_PRECISION( solution, p->v_start, p->v_end, l, threading );

          printf0("g (proc=%d) ---> approx rel error BJ = %f\n", g.my_rank, tmpx1/tmpx2);

          free(solution);
          free(rhs);
          free(x);
          free(diff_sol);
        }
      }

      //MPI_Finalize();
      //exit(0);
      */
      // --------------------------------------------------------------------------------
      //p->eval_operator( output, input, p->op, l, threading );
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
