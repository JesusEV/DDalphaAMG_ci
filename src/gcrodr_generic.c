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

#include "main.h"


// build A and B for the generalized eigenvalue problem within FL-GCRO-DR
void gev_buildAB_PRECISION(complex_PRECISION **A, complex_PRECISION **B, complex_PRECISION **G,
                           vector_PRECISION *W, vector_PRECISION *Z, int mk){

  // TODO

}


int fgmresx_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* This is a trimmed version of FGMRES, as a helper to flgcrodr_PRECISION(...)
* This function assumes:
*	-- p->r already contains the computed residual
*********************************************************************************/  

  // start and end indices for vector functions depending on thread
  int start;
  int end;

  int j=-1, finish=0, iter=0, il;

  PRECISION norm_r0=1, gamma_jp1=1;

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
  
  norm_r0 = creal(p->gamma[0]);

  vector_PRECISION_real_scale( p->V[0], p->r, 1/p->gamma[0], start, end, l ); // v_0 = r / gamma_0
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
  if ( l->level == 0 && l->depth > 0 ) {
    arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, 0, p->preconditioner, p, l, threading );
  }
#endif   
    
  for( il=0; il<p->restart_length && finish==0; il++) {
    j = il; iter++;
      
    // one step of Arnoldi
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      if ( !arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j+1, p->preconditioner, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+2, j+1 );
        break;
      }
    } else {
      if ( !arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
        break;
      }
    }
#else
    if ( !arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p, l, threading ) ) {
      printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
      break;
    }
#endif
      
    if ( cabs( p->H[j][j+1] ) > p->tol/10 ) {
      qr_update_PRECISION( p->H, p->s, p->c, p->gamma, j, l, threading );
      gamma_jp1 = cabs( p->gamma[j+1] );
/*
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
      if ( iter%10 == 0 || p->preconditioner != NULL || l->depth > 0 ) {
        START_MASTER(threading)
        if ( p->print && g.print > 0 )
          printf0("| approx. rel. res. after  %-6d iterations: %e |\n", iter, gamma_jp1/norm_r0 );
        END_MASTER(threading)
      }
#endif
*/

      if( gamma_jp1/norm_r0 < p->tol || gamma_jp1/norm_r0 > 1E+5 ) { // if satisfied ... stop
        
        // FIXME : for now, don't set <finish> i.e. do full FGMRES cycles
        //finish = 1;
        
        START_MASTER(threading)
        if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmres_PRECISION, iter = %d, level=%d\n", iter, l->level );
        END_MASTER(threading)
      }
    } else {
      printf0("depth: %d, iter: %d, p->H(%d,%d) = %+lf+%lfi\n", l->depth, iter, j+1, j, CSPLIT( p->H[j][j+1] ) );
      finish = 1;
      break;
    }
  } // end of the (only and) single restart
  
  if ( l->level == 0 ) {
    START_LOCKED_MASTER(threading)
    g.coarse_iter_count += iter;
    END_LOCKED_MASTER(threading)
  }

  return iter;
}


void flgcrodr_PRECISION_struct_init( gmres_PRECISION_struct *p ) {

  fgmres_PRECISION_struct_init( p );

  // TODO : add here extra stuff associated to FL-GCRO-DR

}


void flgcrodr_PRECISION_struct_alloc( int m, int n, long int vl, PRECISION tol, const int type, const int prec_kind,
                                      void (*precond)(), void (*eval_op)(), gmres_PRECISION_struct *p, level_struct *l ) {

  fgmres_PRECISION_struct_alloc( m, n, vl, tol, type, prec_kind, precond, eval_op, p, l );

  // TODO : add here extra stuff associated to FL-GCRO-DR

}


void flgcrodr_PRECISION_struct_free( gmres_PRECISION_struct *p, level_struct *l ) {

  fgmres_PRECISION_struct_free( p, l );

  // TODO : add here extra stuff associated to FL-GCRO-DR

}


int flgcrodr_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ){

  //printf0("FROM WITHIN GCRODR, depth=%d\n", l->depth);

  // TODO : add <extra> profiling (this extra profiling is to be added to FGMRES)

  // start and end indices for vector functions depending on thread
  int start;
  int end;

  //int res;

  complex_PRECISION beta=0;
  //PRECISION b_norm=1;

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  /*
  PRECISION t0=0, t1=0;

  START_LOCKED_MASTER(threading)
  if ( l->depth==0 && ( p->timing || p->print ) ) prof_init( l );

  if ( l->level==0 && g.num_levels > 1 && g.interpolation ) p->tol = g.coarse_tol;
  if ( l->depth > 0 ) p->timing = 1;
  if ( l->depth == 0 ) t0 = MPI_Wtime();
#if defined(TRACK_RES) && !defined(WILSON_BENCHMARK)
  if ( p->print && g.print > 0 ) printf0("+----------------------------------------------------------+\n");
#endif
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  */























  // ------- FL-GCRO-DR -------------------------------------

  // compute initial residual
  if( p->initial_guess_zero ) {
    //res = _NO_RES;
    vector_PRECISION_copy( p->r, p->b, start, end, l );
  } else {
    //res = _RES;
    if ( p->kind == _LEFT && p->preconditioner ) {
      apply_operator_PRECISION( p->Z[0], p->x, p, l, threading );
      p->preconditioner( p->w, NULL, p->Z[0], _NO_RES, l, threading );
    } else {
      apply_operator_PRECISION( p->w, p->x, p, l, threading ); // compute w = D*x
    }
    vector_PRECISION_minus( p->r, p->b, p->w, start, end, l ); // compute r = b - w
  }

  beta = global_norm_PRECISION( p->r, p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)
  START_MASTER(threading)
  // setting the following line for the upcoming call to fgmresx_PRECISION(...)
  p->gamma[0] = beta;
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);

  if (!p->initial_guess_zero) {
    p->gcrodr_PRECISION.b_norm = global_norm_PRECISION( p->b, p->v_start, p->v_end, l, threading );
    printf0("| initial guess relative residual:            %le |\n", creal(beta)/p->gcrodr_PRECISION.b_norm);
  } else {
    p->gcrodr_PRECISION.b_norm = creal(beta);
  }

  printf0("initial relative residual = %f\n", creal(beta)/p->gcrodr_PRECISION.b_norm);

  if ( p->gcrodr_PRECISION.CU_usable==1 ) {
    // TODO : add code here to build C based on known U
  } else if ( p->gcrodr_PRECISION.CU_usable==0 ) {

    // TODO : add code here to build C and U from scractch
    
    // call one cycle of FGMRES
    fgmresx_PRECISION(p, l, threading);

    complex_PRECISION **A, **B;
    gev_buildAB_PRECISION(A, B, p->H, p->V, p->Z, p->restart_length);

  } else{ error0("Invalid value for p->gcrodr_PRECISION.CU_usable \n"); }

  // --------------------------------------------------------


















  /*
  START_LOCKED_MASTER(threading)
  if ( l->depth == 0 ) { t1 = MPI_Wtime(); g.total_time = t1-t0; g.iter_count = fgmres_iter; }
  END_LOCKED_MASTER(threading)

  if ( l->depth == 0 && g.vt.p_end != NULL  ) {
    if ( g.vt.p_end != NULL ) {
      START_LOCKED_MASTER(threading)
      printf0("solve iter: %d\n", fgmres_iter );
      printf0("solve time: %le seconds\n", t1-t0 );
      g.vt.p_end->values[_SLV_TIME] += (t1-t0)/((double)g.vt.average_over);
      g.vt.p_end->values[_SLV_ITER] += fgmres_iter/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_ITER] += (((double)g.coarse_iter_count)/((double)fgmres_iter))/((double)g.vt.average_over);
      g.vt.p_end->values[_CRS_TIME] += g.coarse_time/((double)g.vt.average_over);
    END_LOCKED_MASTER(threading)
    }
  }
  if ( l->depth == 0 && ( p->timing || p->print ) && !(g.vt.p_end != NULL )  ) {
    START_MASTER(threading)
    if ( g.method != 6 ) prof_print( l );
    END_MASTER(threading)
  }
  */

  int fgmres_iter=0;
  fgmres_iter = fgmres_PRECISION( p, l, threading );

  return fgmres_iter;
}
