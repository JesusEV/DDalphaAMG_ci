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

            t0 = MPI_Wtime();
//	    mumps_setup_PRECISION(lx, threading);

               


#define ICNTL(I) icntl[(I) -1]	//macro according to docu //bridges from fortran indices to c

//            CMUMPS_STRUC_C mumps_id;
            int i;

//######### SET UP RHS #############
            int rhs_len = lx->p_PRECISION.v_end-lx->p_PRECISION.v_start;
            for (i = 0; i < rhs_len; i++){	//set the rhs-indices to global values
              *(px->mumps_irhs_loc + i) = g.my_rank * rhs_len + i+1;		//+1 due to fortran indexing
            }
            vector_PRECISION_copy(px->mumps_rhs_loc, px->b, px->v_start, px->v_end, lx ); //save copy since mumps may change content

 //centralized solution, definitely mention it!
            if (g.my_rank == 0){
              g.mumps_id.rhs = px->mumps_SOL;
            }
            mumps_setup_time = MPI_Wtime() - t0;

	/*
            t0 = MPI_Wtime();
            g.mumps_id.job = 4;	//analyze factorize
            cmumps_c(&(g.mumps_id));
            mumps_job4_time = MPI_Wtime() - t0;
	*/
            mumps_job4_time = -1;
            // 3. solve!
            t0 = MPI_Wtime();
            g.mumps_id.job = 3;		// solve
            cmumps_c(&(g.mumps_id));
            mumps_job3_time = MPI_Wtime() - t0;

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

            int send_count = (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start);
            MPI_Scatter(px->mumps_SOL, send_count, MPI_COMPLEX_PRECISION, px->x, send_count, MPI_COMPLEX_PRECISION, 0, MPI_COMM_WORLD);

            vector_PRECISION mumps_eta=NULL;	//will contain the product of A * Sol_dist
            START_MASTER(threading)
              MALLOC(mumps_eta, complex_PRECISION, lx->vector_size);
            END_MASTER(threading)
            memset(mumps_eta, 0, (lx->p_PRECISION.v_end - lx->p_PRECISION.v_start) * sizeof(complex_PRECISION));


            apply_coarse_operator_PRECISION(mumps_eta, px->x, px->op, lx, threading );
            vector_PRECISION_minus(mumps_eta, mumps_eta, px->b, 0, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start), lx);
            PRECISION mumps_norm = global_norm_PRECISION( mumps_eta, 0, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start), lx, threading );
            PRECISION b_norm = global_norm_PRECISION(px->b, 0, (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start), lx, threading );
            mumps_verify_time = MPI_Wtime() - t0;


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
outfile = fopen("timing_3.txt", "a");
//fprintf(outfile, "Msetup: %f, Ma+f: %f, Msolve: %f, Mverify: %f, FGMRES: %f, FGMRES_iter: %i, BLR: %e, rr: %6f, no_processes: %i\n", mumps_setup_time, mumps_job4_time, mumps_job3_time, mumps_verify_time, t1- t0, nr_iters_gmres, g.mumps_id.cntl[6] ,(mumps_norm/b_norm), g.num_processes);
printf0("Msetup: %f, Ma+f: %f, Msolve: %f, Mverify: %f, FGMRES: %f, FGMRES_iter: %i, BLR: %e, rr: %6f, no_processes: %i\n", mumps_setup_time, mumps_job4_time, mumps_job3_time, mumps_verify_time, t1- t0, nr_iters_gmres, g.mumps_id.cntl[6] ,(mumps_norm/b_norm), g.num_processes);
fclose(outfile);
}
MPI_Barrier(MPI_COMM_WORLD);


//            printf0("Ma+f: %f, Msolve: %f, FGMRES: %f, FGMRES_iter %i, BLR: %e, rr: %6f\n", mumps_job4_time, mumps_job3_time, t1- t0, nr_iters_gmres, g.mumps_id.cntl[6] ,(mumps_norm/b_norm));




            MPI_Scatter(px->mumps_SOL, send_count, MPI_COMPLEX_PRECISION, px->x, send_count, MPI_COMPLEX_PRECISION, 0, MPI_COMM_WORLD);	//scatter again to have px->x filled with mumps' solution
            // ---------------------------------------------------------------------------
/*
printf0("stopping....\n");
MPI_Finalize();
exit(0);
*/
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
