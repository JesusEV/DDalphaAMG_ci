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

#ifdef GCRODR


// declarations of aux functions
int  fgmresx_PRECISION( gmres_PRECISION_struct*, level_struct*, struct Thread* );
void gev_buildAB_PRECISION(complex_PRECISION**, complex_PRECISION**, complex_PRECISION**, vector_PRECISION*, vector_PRECISION*,
                           int, gmres_PRECISION_struct*, level_struct*, struct Thread*);
void order_pairs_PRECISION( vector_PRECISION, complex_PRECISION*, int*, int );
void build_CU_PRECISION( complex_PRECISION**, vector_PRECISION*, vector_PRECISION*,
                         gmres_PRECISION_struct*, level_struct*, struct Thread*, int );
int  arnoldix_step_PRECISION( vector_PRECISION *V, vector_PRECISION *Z, vector_PRECISION w,
                              complex_PRECISION **H, complex_PRECISION* buffer, int j, void (*prec)(),
                              gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading );

// ------------------------------------------------------------------------------------------------------------------------


void flgcrodr_PRECISION_struct_init( gmres_PRECISION_struct *p ) {

  fgmres_PRECISION_struct_init( p );

  p->gcrodr_PRECISION.eigslvr.ordr_idxs = NULL;
  p->gcrodr_PRECISION.eigslvr.ordr_keyscpy = NULL;
}


void flgcrodr_PRECISION_struct_alloc( int m, int n, long int vl, PRECISION tol, const int type, const int prec_kind,
                                      void (*precond)(), void (*eval_op)(), gmres_PRECISION_struct *p, level_struct *l ) {

  fgmres_PRECISION_struct_alloc( m, n, vl, tol, type, prec_kind, precond, eval_op, p, l );

  // g_ln is the length m+k of subspaces used in FL-GCRO-DR
  int g_ln = p->restart_length + p->gcrodr_PRECISION.k;
  int i;
  MALLOC( p->gcrodr_PRECISION.Bbuff, complex_PRECISION*, g_ln );

  p->gcrodr_PRECISION.Bbuff[0] = NULL;
  MALLOC( p->gcrodr_PRECISION.Bbuff[0], complex_PRECISION, g_ln*(g_ln+1) );
  for ( i=1; i<g_ln; i++ ) {
    p->gcrodr_PRECISION.Bbuff[i] = p->gcrodr_PRECISION.Bbuff[0] + i*(g_ln+1);
  }

  MALLOC( p->gcrodr_PRECISION.eigslvr.ordr_idxs, int, g_ln );
  MALLOC( p->gcrodr_PRECISION.eigslvr.ordr_keyscpy, complex_PRECISION, g_ln );
}


void flgcrodr_PRECISION_struct_free( gmres_PRECISION_struct *p, level_struct *l ) {

  fgmres_PRECISION_struct_free( p, l );

  // g_ln is the length m+k of subspaces used in FL-GCRO-DR
  int g_ln = p->restart_length + p->gcrodr_PRECISION.k;
  FREE( p->gcrodr_PRECISION.Bbuff[0], complex_PRECISION, g_ln*(g_ln+1) );
  FREE( p->gcrodr_PRECISION.Bbuff, complex_PRECISION*, g_ln );

  FREE( p->gcrodr_PRECISION.eigslvr.ordr_idxs, int, g_ln );
  FREE( p->gcrodr_PRECISION.eigslvr.ordr_keyscpy, complex_PRECISION, g_ln );
}


int flgcrodr_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ){

  // TODO : add <extra> profiling (this extra profiling is to be added to FGMRES)

  // start and end indices for vector functions depending on thread
  int start;
  int end;

  int res;

  complex_PRECISION beta=0;

  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  // compute initial residual
  if( p->initial_guess_zero ) {
    res = _NO_RES;
    vector_PRECISION_copy( p->r, p->b, start, end, l );
  } else {
    res = _RES;
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

  int fgmresx_iter = 0;

  if ( p->gcrodr_PRECISION.CU_usable==1 ) {
    // TODO : add code here to build C based on known U
  } else if ( p->gcrodr_PRECISION.CU_usable==0 ) {
    // call one cycle of FGMRES
    l->dup_H = 1;
    fgmresx_iter = fgmresx_PRECISION(p, l, threading);
    l->dup_H = 0;

    // update the solution p->x
    compute_solution_PRECISION( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
                                p->y, p->gamma, p->H, fgmresx_iter-1, (res==_NO_RES)?0:1, p, l, threading );

    // TODO : add an update for r here ? (in a relatively cheap way)

    // check if this first call to fgmresx_PRECISION was enough
    PRECISION quot = cabs(p->gamma[fgmresx_iter])/creal(beta);
    if( quot < p->tol || quot > 1E+5 ) {
      START_MASTER(threading)
      if ( quot > 1E+5 )
        printf0("Divergence of fgmresx_PRECISION, iter = %d, level=%d\n", fgmresx_iter, l->level );
      END_MASTER(threading)

      return fgmresx_iter;
    }

    if ( p->preconditioner==NULL ) {
    } else {
    }

    if ( p->preconditioner==NULL ) {
      // build the matrices A and B used for generalized-eigensolving
      gev_buildAB_PRECISION( p->gcrodr_PRECISION.gev_A, p->gcrodr_PRECISION.gev_B, p->gcrodr_PRECISION.eigslvr.Hc,
                             p->V, p->V, fgmresx_iter, p, l, threading );
      // build C and U
      build_CU_PRECISION( p->gcrodr_PRECISION.eigslvr.Hc, p->V, p->V, p, l, threading, fgmresx_iter );
    } else {
      // build the matrices A and B used for generalized-eigensolving
      gev_buildAB_PRECISION( p->gcrodr_PRECISION.gev_A, p->gcrodr_PRECISION.gev_B, p->gcrodr_PRECISION.eigslvr.Hc,
                             p->V, p->Z, fgmresx_iter, p, l, threading );
      // build C and U
      build_CU_PRECISION( p->gcrodr_PRECISION.eigslvr.Hc, p->V, p->Z, p, l, threading, fgmresx_iter );
    }

  } else{ error0("Invalid value for p->gcrodr_PRECISION.CU_usable \n"); }

  // TODO : add main loop for FL-GCRO-DR

  // --------------------------------------------------------

  int fgmres_iter=0;
  p->initial_guess_zero = 0;
  fgmres_iter = fgmres_PRECISION( p, l, threading );
  p->initial_guess_zero = 1;

  return fgmres_iter+fgmresx_iter;
}




// ---------- AUXILIARY FUNCTIONS

// build A and B for the generalized eigenvalue problem within FL-GCRO-DR
void gev_buildAB_PRECISION( complex_PRECISION **A, complex_PRECISION **B, complex_PRECISION **G,
                            vector_PRECISION *W, vector_PRECISION *Z, int mk, gmres_PRECISION_struct *p,
                            level_struct *l, struct Thread *threading ){

  int i,j,k;
  complex_PRECISION **Bbuff = p->gcrodr_PRECISION.Bbuff;
  complex_PRECISION tmp[mk+1];

  // -------- building B

  for ( j=0; j<mk; j++ ) {
    process_multi_inner_product_PRECISION( mk+1, Bbuff[j], W, Z[j], p->v_start, p->v_end, l, threading );

    START_MASTER(threading)
    for( i=0; i<(mk+1); i++ )
      tmp[i] = Bbuff[j][i];
    if ( g.num_processes > 1 ) {
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( tmp, Bbuff[j], mk+1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

  }

  // FIXME : improve the following matrix-matrix multiplication by using more threads than <master>
  START_MASTER(threading)
  for ( j=0; j<mk; j++ ) {
    for ( i=0; i<mk; i++ ) {
      B[j][i] = 0.0;
      for ( k=0; k<mk+1; k++ ) {
        B[j][i] += conj_PRECISION(G[i][k])*(Bbuff[j][k]);
      }
    }
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // -------- building A

  // FIXME : improve the following matrix-matrix multiplication by using more threads than <master>
  START_MASTER(threading)
  for ( j=0; j<mk; j++ ) {
    for ( i=0; i<mk; i++ ) {
      A[j][i] = 0.0;
      for ( k=0; k<mk+1; k++ ) {
        A[j][i] += conj_PRECISION(G[i][k])*(G[j][k]);
      }
    }
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
}


int fgmresx_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* This is a trimmed version of FGMRES, as a helper to flgcrodr_PRECISION(...)
* This function assumes:
*	-- p->r already contains the computed residual
*********************************************************************************/  

  // TODO : what happens if a <breakdown> happens, and we haven't reached k iters ??

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
    arnoldix_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, 0, p->preconditioner, p, l, threading );
  }
#endif   
    
  for( il=0; il<p->restart_length && finish==0; il++) {
    j = il; iter++;
      
    // one step of Arnoldi
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    if ( l->level == 0 && l->depth > 0 ) {
      if ( !arnoldix_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j+1, p->preconditioner, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+2, j+1 );
        break;
      }
    } else {
      if ( !arnoldix_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p, l, threading ) ) {
        printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
        break;
      }
    }
#else
    if ( !arnoldix_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j, p->preconditioner, p, l, threading ) ) {
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

      // do at least GCRODR.k iters, to be able to construct the recycling subspace
      if( (gamma_jp1/norm_r0 < p->tol || gamma_jp1/norm_r0 > 1E+5) && ((j+1) >= p->gcrodr_PRECISION.k) ) { // if satisfied ... stop
        
        finish = 1;
        
        START_MASTER(threading)
        if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmresx_PRECISION, iter = %d, level=%d\n", iter, l->level );
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


// order the eigenpairs returned by the generalized eigenvalue solver from LAPACK
void order_pairs_PRECISION( vector_PRECISION keys, complex_PRECISION *keys_cpy, int *output, int n ){
  int i,j;
  int buff1;

  complex_PRECISION buff2;

  memcpy(keys_cpy, keys, sizeof(complex_PRECISION)*n);

  for ( i=0; i<n; i++ ) {
    output[i] = i;
  }

  for ( i=0; i<n; i++ ) {
    for ( j=i+1; j<n; j++ ) {
      if ( cabs(keys_cpy[i]) > cabs(keys_cpy[j]) ) {
        buff2 =  keys_cpy[i];
        keys_cpy[i] = keys_cpy[j];
        keys_cpy[j] = buff2;

        buff1 =  output[i];
        output[i] = output[j];
        output[j] = buff1;
      }
    }
  }
}


void build_CU_PRECISION( complex_PRECISION **G, vector_PRECISION *W, vector_PRECISION *Z,
                         gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading, int m ){

  // --------- eigensolving first

  // FIXME : improve the following eigensolve by using more threads than <master>
  START_MASTER(threading)
  // calling LAPACK's generalized eigenvalue solver through LAPACKE
  p->gcrodr_PRECISION.eigslvr.N = m;
  gen_eigslvr_PRECISION( &(p->gcrodr_PRECISION.eigslvr) );
  // p->gcrodr_PRECISION.eigslvr.ordr_idxs contains the indices to access w and vr in ascending magnitude
  order_pairs_PRECISION( p->gcrodr_PRECISION.eigslvr.w, p->gcrodr_PRECISION.eigslvr.ordr_keyscpy,
                         p->gcrodr_PRECISION.eigslvr.ordr_idxs, m );
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // ---------------- then, computing C and U

  int i, j, kl, start, end;
  compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);

  int g_ln = p->restart_length + p->gcrodr_PRECISION.k;

  vector_PRECISION *Yk = p->gcrodr_PRECISION.Yk;
  vector_PRECISION *Ck = p->gcrodr_PRECISION.C;
  vector_PRECISION *Uk = p->gcrodr_PRECISION.U;
  //vector_PRECISION *Z  = (p->preconditioner==NULL) ? p->V : p->Z;
  vector_PRECISION vr  = p->gcrodr_PRECISION.eigslvr.vr;

  int *idxs = p->gcrodr_PRECISION.eigslvr.ordr_idxs;
  int k = p->gcrodr_PRECISION.k;

  complex_PRECISION **QR = p->gcrodr_PRECISION.QR;
  complex_PRECISION **Q = p->gcrodr_PRECISION.Q;
  complex_PRECISION **R = p->gcrodr_PRECISION.R;
  complex_PRECISION **Rinv = p->gcrodr_PRECISION.eigslvr.qr_Rinv;

  // for each new eigensolution, we have a new mapping for Pk
  vector_PRECISION *Pk = p->gcrodr_PRECISION.Pk;
  for ( i=0; i<k; i++ ) {
    Pk[i] = vr + idxs[i] * g_ln;
  }

  // compute Yk
  for ( i=0; i<k; i++ ) {
    // set all vectors in Yk to zero, to accumulate
    vector_PRECISION_define( Yk[i], 0, start, end, l );
    // and then, multi saxpy to obtain Yk
    vector_PRECISION_multi_saxpy( Yk[i], Z, Pk[i], 1, m, start, end, l );
  }

  // build the matrix for QR
  // FIXME : improve the following matrix-matrix multiplication by using more threads than <master>
  START_MASTER(threading)
  for ( j=0; j<k; j++ ) {
    // set all column j of QR to zero
    memset( QR[j], 0.0, sizeof(complex_PRECISION)*(m+1) );
    // and then do accumulations over that column
    for ( i=0; i<m; i++ ) {
      for ( kl=0; kl<(m+1); kl++ ) {
        QR[j][kl] += G[i][kl] * Pk[j][i];
      }
    }
  }
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // QR decomposition of G*Pk
  // FIXME : improve the following matrix-matrix multiplication by using more threads than <master>
  START_MASTER(threading)

  // QR
  p->gcrodr_PRECISION.eigslvr.qr_m = m+1;
  p->gcrodr_PRECISION.eigslvr.qr_n = k;
  qr_PRECISION( &(p->gcrodr_PRECISION.eigslvr) );

  p->gcrodr_PRECISION.eigslvr.qr_k = k;

  // compute R^{-1}
  for ( j=0; j<k; j++ ) {
    for ( i=0; i<k; i++ ) {
      R[j][i] = QR[j][i];
    }
  }
  inv_tri_PRECISION( &(p->gcrodr_PRECISION.eigslvr) );

  // compute Q
  q_from_qr_PRECISION( &(p->gcrodr_PRECISION.eigslvr) );

  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)

  // compute Ck
  for ( i=0; i<k; i++ ) {
    // set all vectors in Yk to zero, to accumulate
    vector_PRECISION_define( Ck[i], 0, start, end, l );
    // and then, multi saxpy to obtain Yk
    vector_PRECISION_multi_saxpy( Ck[i], W, Q[i], 1, m+1, start, end, l );
  }

  // compute Uk
  for ( i=0; i<k; i++ ) {
    // set all vectors in Yk to zero, to accumulate
    vector_PRECISION_define( Uk[i], 0, start, end, l );
    // and then, multi saxpy to obtain Yk
    vector_PRECISION_multi_saxpy( Uk[i], Yk, Rinv[i], 1, k, start, end, l );
  }
}


int arnoldix_step_PRECISION( vector_PRECISION *V, vector_PRECISION *Z, vector_PRECISION w,
                             complex_PRECISION **H, complex_PRECISION* buffer, int j, void (*prec)(),
                             gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

  // TODO : extend this to include orthonormalization against p->gcrodr_PRECISION.C

  int return_val;
  return_val = arnoldi_step_PRECISION( V, Z, w, H, buffer, j, prec, p, l, threading );
  return return_val;

}

#endif
