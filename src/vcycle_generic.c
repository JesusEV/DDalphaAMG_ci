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
#include "vcycle_PRECISION.h"

void smoother_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                         int n, const int res, level_struct *l, struct Thread *threading ) {
  
  ASSERT( phi != eta );

  START_MASTER(threading);
  PROF_PRECISION_START( _SM );
  END_MASTER(threading);
  
  if ( g.method == 1 ) {
    additive_schwarz_PRECISION( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else if ( g.method == 2 ) {
    red_black_schwarz_PRECISION( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else if ( g.method == 3 ) {
    sixteen_color_schwarz_PRECISION( phi, Dphi, eta, n, res, &(l->s_PRECISION), l, threading );
  } else {
    int start = threading->start_index[l->depth];
    int end   = threading->end_index[l->depth];
    START_LOCKED_MASTER(threading)
    l->sp_PRECISION.initial_guess_zero = res;
    l->sp_PRECISION.num_restart = n;
    END_LOCKED_MASTER(threading)
    if ( g.method == 4 || g.method == 6 ) {
      if ( g.odd_even ) {
        if ( res == _RES ) {
          apply_operator_PRECISION( l->sp_PRECISION.x, phi, &(l->p_PRECISION), l, threading );
          vector_PRECISION_minus( l->sp_PRECISION.x, eta, l->sp_PRECISION.x, start, end, l );
        }
        block_to_oddeven_PRECISION( l->sp_PRECISION.b, res==_RES?l->sp_PRECISION.x:eta, l, threading );
        START_LOCKED_MASTER(threading)
        l->sp_PRECISION.initial_guess_zero = _NO_RES;
        END_LOCKED_MASTER(threading)
        if ( g.method == 6 ) {
          if ( l->depth == 0 ) g5D_solve_oddeven_PRECISION( &(l->sp_PRECISION), &(l->oe_op_PRECISION), l, threading );
          else g5D_coarse_solve_odd_even_PRECISION( &(l->sp_PRECISION), &(l->oe_op_PRECISION), l, threading );
        } else {
          if ( l->depth == 0 ) solve_oddeven_PRECISION( &(l->sp_PRECISION), &(l->oe_op_PRECISION), l, threading );
          else coarse_solve_odd_even_PRECISION( &(l->sp_PRECISION), &(l->oe_op_PRECISION), l, threading );
        }
        if ( res == _NO_RES ) {
          oddeven_to_block_PRECISION( phi, l->sp_PRECISION.x, l, threading );
        } else {
          oddeven_to_block_PRECISION( l->sp_PRECISION.b, l->sp_PRECISION.x, l, threading );
          vector_PRECISION_plus( phi, phi, l->sp_PRECISION.b, start, end, l );
        }
      } else {
        START_LOCKED_MASTER(threading)
        l->sp_PRECISION.x = phi; l->sp_PRECISION.b = eta;
        END_LOCKED_MASTER(threading)
        fgmres_PRECISION( &(l->sp_PRECISION), l, threading );
      }
    } else if ( g.method == 5 ) {
      vector_PRECISION_copy( l->sp_PRECISION.b, eta, start, end, l );
      bicgstab_PRECISION( &(l->sp_PRECISION), l, threading );
      vector_PRECISION_copy( phi, l->sp_PRECISION.x, start, end, l );
    }
    ASSERT( Dphi == NULL );
  }
  
  START_MASTER(threading);
  PROF_PRECISION_STOP( _SM, n );
  END_MASTER(threading);
}


void vcycle_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                       int res, level_struct *l, struct Thread *threading ) {

  if ( g.interpolation && l->level>0 ) {
    for ( int i=0; i<l->n_cy; i++ ) {
      if ( i==0 && res == _NO_RES ) {
        restrict_PRECISION( l->next_level->p_PRECISION.b, eta, l, threading );
      } else {
        int start = threading->start_index[l->depth];
        int end   = threading->end_index[l->depth];
        apply_operator_PRECISION( l->vbuf_PRECISION[2], phi, &(l->p_PRECISION), l, threading );
        vector_PRECISION_minus( l->vbuf_PRECISION[3], eta, l->vbuf_PRECISION[2], start, end, l );
        restrict_PRECISION( l->next_level->p_PRECISION.b, l->vbuf_PRECISION[3], l, threading );
      }
      if ( !l->next_level->idle ) {
        START_MASTER(threading)
        if ( l->depth == 0 )
          g.coarse_time -= MPI_Wtime();
        END_MASTER(threading)
        if ( l->level > 1 ) {
          if ( g.kcycle )
            fgmres_PRECISION( &(l->next_level->p_PRECISION), l->next_level, threading );
          else
            vcycle_PRECISION( l->next_level->p_PRECISION.x, NULL, l->next_level->p_PRECISION.b, _NO_RES, l->next_level, threading );
        } else {
          if ( g.odd_even ) {
            if ( g.method == 6 ) {
              g5D_coarse_solve_odd_even_PRECISION( &(l->next_level->p_PRECISION), &(l->next_level->oe_op_PRECISION), l->next_level, threading );
            } else {

              START_MASTER(threading)
              g.coarsest_time -= MPI_Wtime();
              END_MASTER(threading)

              coarse_solve_odd_even_PRECISION( &(l->next_level->p_PRECISION), &(l->next_level->oe_op_PRECISION), l->next_level, threading );

              START_MASTER(threading)
              g.coarsest_time += MPI_Wtime();
              END_MASTER(threading)

            }
          } else {

            // ---------------------------------------------------------------------------
            gmres_PRECISION_struct* px = &(l->next_level->p_PRECISION);
            level_struct* lx = l->next_level;

            double t0,t1;	//timing
            double mumps_setup_time, mumps_job4_time, mumps_job3_time, mumps_verify_time;




            int num_eig_vect = lx->num_parent_eig_vect, 
            vector_size = lx->num_lattice_site_var,
            clover_size = (2*num_eig_vect*num_eig_vect+num_eig_vect);

            int start, end;
            compute_core_start_end_custom(0, lx->num_inner_lattice_sites, &start, &end, lx, threading, 1);
            // 1. prepare sparse matrix for MUMPS (take into account : matrix-vector mult receives : px->op)


            t0 = MPI_Wtime();
	    mumps_setup_PRECISION(lx, threading);

            // 2. analyze+factorize
//############################## MUMPS RELATED STUFF ###############################################

#define ICNTL(I) icntl[(I) -1]	//macro according to docu //bridges from fortran indices to c
//#define CNTL(I) cntl[(I) -1]	//same macro for cntl <-- does not work, messes up the line above.

//            CMUMPS_STRUC_C mumps_id;
            int i;
            int chunklen = SQUARE(lx->num_lattice_site_var) *lx->num_inner_lattice_sites *9;

            int mumps_n = lx->num_lattice_site_var * lx->num_inner_lattice_sites * g.num_processes;	//order of Matrix
            int nnz = chunklen * g.num_processes;	//number of nonzero elements
            int nnz_loc = chunklen;


//######### SET UP RHS #############

//vector_PRECISION_define_random(px->b, px->v_start, px->v_end, lx );


            int rhs_len = lx->p_PRECISION.v_end-lx->p_PRECISION.v_start;  
            for (i = 0; i < rhs_len; i++){	//set the rhs-indices to global values
              *(px->mumps_irhs_loc + i) = g.my_rank * rhs_len + i+1;		//+1 due to fortran indexing
              *(px->mumps_rhs_loc + i) = px->b[i];				//save copy since mumps may change content
            }
 //centralized solution, definitely mention it!
            complex_PRECISION* SOL;
            if (g.my_rank == 0){
              START_MASTER(threading)
                MALLOC(SOL, complex_PRECISION, mumps_n);
              END_MASTER(threading)
              memset(SOL, 0, mumps_n * sizeof(complex_PRECISION));
              g.mumps_id.rhs = SOL;
            }
            mumps_setup_time = MPI_Wtime() - t0;

	
            t0 = MPI_Wtime();
            g.mumps_id.job = 4;	//analyze factorize
            cmumps_c(&(g.mumps_id));
            mumps_job4_time = MPI_Wtime() - t0;


            // 3. solve!
            t0 = MPI_Wtime();
            g.mumps_id.job = 3;		// solve
            cmumps_c(&(g.mumps_id));
            mumps_job3_time = MPI_Wtime() - t0;



//            MPI_Finalize();
//            exit(0);
           
//####################### CHECK RELATIVE RESIDUAL ##################################
/*
0. Scatter SOL to SOL_dist
1. find A * SOL
2. compute vector_PRECISION_minus A*SOL - eta
3. compute || A*SOL - eta || and || eta ||
4. divide both for relative residual
*/
            t0 = MPI_Wtime();
            int nx = (SQUARE(lx->num_lattice_site_var)*9) * lx->num_inner_lattice_sites;

//	###################### 0. ######################
            vector_PRECISION SOL_dist=NULL;	//will contain distributed Solution
            START_MASTER(threading)
              MALLOC(SOL_dist, complex_PRECISION, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start));
            END_MASTER(threading)
            memset(SOL_dist, 0, (lx->p_PRECISION.v_end - lx->p_PRECISION.v_start) * sizeof(complex_PRECISION));
//			scatter SOL to SOL_dist on each process
            int send_count = (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start);
//MPI_Scatter(void* send_data, int send_count, MPI_Datatype send_datatype, void* recv_data, int recv_count, MPI_Datatype recv_datatype, int root, MPI_Comm communicator)
            MPI_Scatter(SOL, send_count, MPI_COMPLEX_PRECISION, SOL_dist, send_count, MPI_COMPLEX_PRECISION, 0, MPI_COMM_WORLD);

//	###################### 1. ######################
//			allocate and set new "eta" vector
            vector_PRECISION mumps_eta=NULL;	//will contain the product of A * Sol_dist

            START_MASTER(threading)
              MALLOC(mumps_eta, complex_PRECISION, lx->vector_size);
            END_MASTER(threading)
            memset(mumps_eta, 0, (lx->p_PRECISION.v_end - lx->p_PRECISION.v_start) * sizeof(complex_PRECISION));

            int nxy = (SQUARE(lx->num_lattice_site_var)*9) * lx->num_inner_lattice_sites;
//spmv_PRECISION(mumps_eta, SOL_dist, lx->p_PRECISION.mumps_vals, lx->p_PRECISION.mumps_Is, lx->p_PRECISION.mumps_Js, nxy, &(lx->p_PRECISION), lx, threading );

            apply_coarse_operator_PRECISION(mumps_eta, SOL_dist, px->op, lx, threading );
//            apply_operator_PRECISION(mumps_eta, SOL_dist, px, lx, threading );

//	###################### 2. ######################
//void vector_PRECISION_minus( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, int start, int end, level_struct *l ); // z := x - y
//            vector_PRECISION_minus(mumps_eta, mumps_eta, px->b, 0, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start), lx);
            vector_PRECISION_minus(mumps_eta, mumps_eta, px->b, 0, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start), lx);

//	###################### 3. ######################
//PRECISION global_norm_PRECISION( vector_PRECISION phi, int start, int end, level_struct *l, struct Thread *threading );
//  PRECISION mumps_norm;
//            PRECISION mumps_norm = global_norm_PRECISION( mumps_eta, 0, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start), lx, threading );
            PRECISION mumps_norm = global_norm_PRECISION( mumps_eta, 0, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start), lx, threading );
//  PRECISION DD_norm;
            PRECISION b_norm = global_norm_PRECISION(px->b, 0, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start), lx, threading );

//	###################### 4. ######################
  //PRECISION rr;
            mumps_verify_time = MPI_Wtime() - t0;
  //          printf0("rr: %6f,\tmumps_norm:%6f,\tb_norm:%6f\n", (mumps_norm/b_norm), mumps_norm, b_norm);

// finish mumps
            t0 = MPI_Wtime();
//            g.mumps_id.job = JOB_END;
//            cmumps_c(&(g.mumps_id));
            mumps_setup_time += MPI_Wtime() - t0;


            t0 = MPI_Wtime();
            int nr_iters_gmres = fgmres_PRECISION(px, lx, threading );
				//solution will be in px->x
            t1 = MPI_Wtime();

if (g.my_rank == 0){
FILE *outfile;
outfile = fopen("timing_2.txt", "a");
//fprintf(outfile, "Msetup: %f, Ma+f: %f, Msolve: %f, Mverify: %f, FGMRES: %f, FGMRES_iter: %i, BLR: %e, rr: %6f, no_processes: %i\n", mumps_setup_time, mumps_job4_time, mumps_job3_time, mumps_verify_time, t1- t0, nr_iters_gmres, g.mumps_id.cntl[6] ,(mumps_norm/b_norm), g.num_processes);
printf0("Msetup: %f, Ma+f: %f, Msolve: %f, Mverify: %f, FGMRES: %f, FGMRES_iter: %i, BLR: %e, rr: %6f, no_processes: %i\n", mumps_setup_time, mumps_job4_time, mumps_job3_time, mumps_verify_time, t1- t0, nr_iters_gmres, g.mumps_id.cntl[6] ,(mumps_norm/b_norm), g.num_processes);
fclose(outfile);
}
MPI_Barrier(MPI_COMM_WORLD);


//            printf0("Ma+f: %f, Msolve: %f, FGMRES: %f, FGMRES_iter %i, BLR: %e, rr: %6f\n", mumps_job4_time, mumps_job3_time, t1- t0, nr_iters_gmres, g.mumps_id.cntl[6] ,(mumps_norm/b_norm));

            printf0("stopping ...\n");
            MPI_Finalize();
            exit(0);

            // ---------------------------------------------------------------------------

          }
        }
        START_MASTER(threading)
        if ( l->depth == 0 )
          g.coarse_time += MPI_Wtime();
        END_MASTER(threading)
      }
      if( i == 0 && res == _NO_RES )
        interpolate3_PRECISION( phi, l->next_level->p_PRECISION.x, l, threading );
      else
        interpolate_PRECISION( phi, l->next_level->p_PRECISION.x, l, threading );
      smoother_PRECISION( phi, Dphi, eta, l->post_smooth_iter, _RES, l, threading );
      res = _RES;
    }
  } else {
    smoother_PRECISION( phi, Dphi, eta, (l->depth==0)?l->n_cy:l->post_smooth_iter, res, l, threading );
  }
}
