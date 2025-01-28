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

extern global_struct g;
#ifdef HAVE_HDF5
Hdf5_fileinfo h5info;
#endif
extern struct common_thread_data *commonthreaddata;
extern struct Thread *no_threading;

int main( int argc, char **argv ) {
    
#ifdef HAVE_HDF5
  h5info.filename=NULL;
  h5info.file_id=-1; 
  h5info.rootgroup_id=-1; 
  h5info.configgroup_id=-1;
  h5info.eigenmodegroup_id=-1;
  h5info.thiseigenmodegroup_id=-1;
  h5info.isOpen=0;
  h5info.mode=-1;
#endif
  level_struct l;
  config_double hopp = NULL;
  
  MPI_Init( &argc, &argv );
  
  predefine_rank( MPI_COMM_WORLD );
  if ( g.my_rank == 0 ) {
    printf("\n\n+----------------------------------------------------------+\n");
    printf("| The DDalphaAMG solver library.                           |\n");
    printf("| Copyright (C) 2016, Matthias Rottmann, Artur Strebel,    |\n");
    printf("|       Simon Heybrock, Simone Bacchio, Bjoern Leder.      |\n");
    printf("|                                                          |\n");
    printf("| This program comes with ABSOLUTELY NO WARRANTY.          |\n");
    printf("+----------------------------------------------------------+\n\n");
  }
  
  method_init( &argc, &argv, &l );
  
  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);
  
  MALLOC( hopp, complex_double, 3*l.inner_vector_size );

  if(g.in_format == _LIME)
    lime_read_conf( (double*)(hopp), g.in, &(g.plaq_hopp) );
  else 
    read_conf( (double*)(hopp), g.in, &(g.plaq_hopp), &l );

  // store configuration, compute clover term
  dirac_setup( hopp, &l );
  FREE( hopp, complex_double, 3*l.inner_vector_size );

  commonthreaddata = (struct common_thread_data *)malloc(sizeof(struct common_thread_data));
  init_common_thread_data(commonthreaddata);
  
  THREADED(g.num_openmp_processes)
  {
    g.on_solve = 0;
    struct Thread threading;
    setup_threading(&threading, commonthreaddata, &l);
    setup_no_threading(no_threading, &l);

    //double t0x=0, t1x=0, elap_time=0;

    //t0x = MPI_Wtime();


    // setup up initial MG hierarchy
    method_setup( NULL, &l, &threading );

    //t1x = MPI_Wtime();
    //elap_time = t1x-t0x;
    //if (g.my_rank==0) printf("elapsed time (init setup phase): %-8.4lf seconds\n", elap_time);

    //t0x = MPI_Wtime();

#ifdef MUMPS_ADDS
    {
      level_struct *lx = &l;
      int i;
      for (i = 1; i<g.num_levels; i++) {
        lx = lx->next_level;
      }
      if (!lx->idle){

        struct Thread* threadx = &threading;

        SYNC_CORES(threadx)
        START_MASTER(threadx)
	printf0("call to mumps_setup from main.c\n");
	END_MASTER(threadx)
        mumps_setup_float(lx, &threading);        //setup vals, Is, Js
        START_MASTER(threadx)
        printf0("mumps_setup done in main.c\n");
	END_MASTER(threadx)
        SYNC_CORES(threadx)

        double t0,t1;
        START_MASTER(threadx)
        t0 = MPI_Wtime();

#ifndef DenseDirectSolves
        printf0("starting analyze from main.c\n");
        END_MASTER(threadx)
        SYNC_CORES(threadx)
    
        START_MASTER(threadx)
        //g.mumps_id.job = 4; //analyze and factorize
        g.mumps_id.job = 1; //analyze

        cmumps_c(&(g.mumps_id));
        printf0("analyze done, starting factorize singlethreaded from main.c\n");
  
  
        g.mumps_id.job = 2; //factorize
        cmumps_c(&(g.mumps_id));
#else
        printf0("starting inverting using scalapack from main.c\n");
//Testing from here

	vector_float rand, xs, xg;
	MALLOC( rand, complex_float, 2*lx->vector_size);
	MALLOC( xs, complex_float, 2*lx->vector_size);
	MALLOC( xg, complex_float, 2*lx->vector_size);
	memset( rand, 0, 2*lx->vector_size * sizeof(complex_float));
	memset( xs, 0, 2*lx->vector_size * sizeof(complex_float));
	memset( xg, 0, 2*lx->vector_size * sizeof(complex_float));

	float r1, r2, r3; //norms 
	int start = 0, end = 2*lx->vector_size;
	vector_float_define_random( rand, start, end, lx );
	vector_float_copy( lx->p_float.b, rand, start, end, lx );	 // r = eta from start to end on level l

//printing matrix
	printf0("D = \n");
	int of = 0; //55
	int size = 3;
	for (int i=0; i < lx->num_inner_lattice_sites * lx->num_lattice_site_var; i = i + size){
	    for (int j=0; j < lx->num_inner_lattice_sites * lx->num_lattice_site_var; j = j + size){
		if (fabs((lx->p_float.dense_vals[(i + of) + lx->num_inner_lattice_sites *
			    lx->num_lattice_site_var * (j + of)])) > 1e-7) printf0("X");
		else printf0(" ");
//		printf0("%e\t", (creal(lx->p_float.dense_vals[i * lx->num_inner_lattice_sites * lx->num_lattice_site_var + j])));
	    }
	    printf0("\n");
	}



	//x = A *b via scalapack currently in "invert_coarsest_matrix" is apply coarsest matrix 
	printf0("applying coarsest matrix via scalapack from main.c\n");
	invert_coarsest_matrix_scalap_float( lx, lx->p_float.dense_vals, lx->p_float.desc_dense_vals, lx->p_float.b, lx->p_float.desc_rhs, g.mumps_id.n, &threading);	
	vector_float_copy( xs, lx->p_float.b, start, end, lx );	 // r = eta from start to end on level l

	//applying matrix via DDalphaAMG
	apply_coarse_operator_float( xg, lx->p_float.b, lx->p_float.op, lx, &threading );

	//checking norm of diff:
	vector_float_minus( rand, xg, xs, start, end, lx );	
	r1 = global_norm_float( rand, start, end, lx, &threading );
	r1 = r1/global_norm_float( xg, start, end, lx, &threading );

	printf0("\n\n\nrelative res. |(Ax)_s - (Ax)_d)|/|(Ax)_d| = %f\n\n\n\n", r1);


//checking elementwise
	printf0("xs\t\t\txg\n");
	for (int i = 0; i < lx->num_inner_lattice_sites * lx->num_lattice_site_var; i = i+28){
	    printf0("%+5.1f%+5.1fi,\t%+5.1f%+5.1fi,\n", CSPLIT(xs[i]), CSPLIT(xg[i]));
	}
	exit(0);
//inverting using scalapack uses p_float.b as RHS and as solution afterwards
	invert_coarsest_matrix_scalap_float( lx, lx->p_float.dense_vals, lx->p_float.desc_dense_vals, lx->p_float.b, lx->p_float.desc_rhs, g.mumps_id.n, &threading);
	
	//compute A * sol_from_scalapack = xs
	printf0("applying coarse operator, in main.c\n");
	//apply_operator_float(xs, lx->p_float.b, lx->p_float.eval_operator, lx, &threading );
	apply_coarse_operator_float( xs, lx->p_float.b, lx->p_float.op, lx, &threading );


	//compute norms of solutions ||rand - xs|| / ||rand||
	printf0("computing relative residual, in main.c\n");

	r2 = global_norm_float( xs, start, end, lx, &threading );
	r3 = global_norm_float( rand, start, end, lx, &threading );
	printf0("\n\n\nrelative res. |A ( A^-1 b)| = %f, \t |b| = %f\n\n\n\n", r2, r3);
	
	vector_float_minus( xs, rand, xs, start, end, lx );	
	r2 = global_norm_float( xs, start, end, lx, &threading );
	
	printf0("\n\n\nrelative res. |A ( A^-1 b) - b| = %f\n\n\n\n", r2);
	
	r2  = r2 / global_norm_float( rand, start, end, lx, &threading );

	printf0("\n\n\nrelative res. |A ( A^-1 b) - b| / |b| = %f\n\n\n\n", r2);

	exit(0);

	//copy scalapack solution to xs
	printf0("copying solution from scalapack to xs, in main.c\n");
	vector_float_copy( xs, lx->p_float.b, start, end, lx );	 // r = eta from start to end on level l
	
	//reset p_float.b to rand
	printf0("resetting b to rand, in main.c\n");
	vector_float_copy( lx->p_float.b, rand, start, end, lx );	 // r = eta from start to end on level l

	//compute solution using fgmres
	printf0("inverting using fgmres, in main.c\n");
	int its = fgmres_float( &(lx->p_float), lx,  &threading );

	//copy fgmres solution to xg
	printf0("copying solution from fgmres to xg, in main.c\n");
	vector_float_copy( xg, lx->p_float.x, start, end, lx );	 // r = eta from start to end on level l

	//compute norms of solutions ||xg - xs|| / ||xg||
	printf0("computing relative residual, in main.c\n");
	vector_float_minus( xs, xg, xs, start, end, lx );	
	r1 = global_norm_float( xs, start, end, lx, &threading );
	r1  = r1 / global_norm_float( xg, start, end, lx, &threading );

	printf0("\n\n\nrelative res. |xg - xs| / |xg| = %f\n\n\n\n", r1);
	//check if p_float.b is changes ||rand - b|| / ||rand||


//	float global_norm_PRECISION( vector_PRECISION phi, int start, int end, level_struct *l, struct Thread *threading );

//Testing to here
	invert_coarsest_matrix_scalap_float( lx, lx->p_float.dense_vals, lx->p_float.desc_dense_vals, lx->p_float.b, lx->p_float.desc_rhs, g.mumps_id.n, &threading);
        printf0("inverting using scalapack done\n");

	exit(0);
#endif

        t1 = MPI_Wtime();
	g.mumps_fact_time += t1 - t0;
#ifndef DenseDirectSolves
	if (g.my_rank == 0) printf("MUMPS analyze and factorize time (seconds) : %f \t in main.c\n",t1-t0);
#else
	if (g.my_rank == 0) printf("Invert using scalapack time (seconds) : %f \t in main.c\n",t1-t0);
#endif
        printf0("factorize done in main.c\n");
        END_MASTER(threadx)
        SYNC_CORES(threadx)
      }
    }
#endif

#if defined(POLYPREC) || defined(GCRODR)
    {
      level_struct *lx = &l;
      while (1) {
        if ( lx->level==0 ) {
          if ( g.mixed_precision==0 ) {
#ifdef GCRODR
            lx->p_double.gcrodr_double.k = g.gcrodr_k_setup;
#endif
#ifdef POLYPREC
            lx->p_float.polyprec_float.d_poly = g.polyprec_d_setup;
#endif
          }
          else {
#ifdef GCRODR
            lx->p_float.gcrodr_float.k = g.gcrodr_k_setup;
#endif
#ifdef POLYPREC
            lx->p_float.polyprec_float.d_poly = g.polyprec_d_setup;
#endif
          }
          break;
        }
        else { lx = lx->next_level; }
      }
    }
#endif

    // iterative phase
    method_update( l.setup_iter, &l, &threading );

    //t1x = MPI_Wtime();
    //elap_time = t1x-t0x;
    //if (g.my_rank==0) printf("elapsed time (iterative setup phase): %-8.4lf seconds\n", elap_time);

#if defined(POLYPREC) || defined(GCRODR)
    {
      level_struct *lx = &l;
      while (1) {
        if ( lx->level==0 ) {
          if ( g.mixed_precision==0 ) {
#ifdef POLYPREC
            lx->p_double.polyprec_double.d_poly = g.polyprec_d_solve;
#endif
#ifdef GCRODR
            lx->p_double.gcrodr_double.k = g.gcrodr_k_solve;
#endif
          }
          else {
#ifdef POLYPREC
            lx->p_float.polyprec_float.d_poly = g.polyprec_d_solve;
#endif
#ifdef GCRODR
            lx->p_float.gcrodr_float.k = g.gcrodr_k_solve;
#endif
          }
          break;
        }
        else { lx = lx->next_level; }
      }
    }
#endif

    g.on_solve = 1;


    solve_driver( &l, &threading );
  }
  
  finalize_common_thread_data(commonthreaddata);
  finalize_no_threading(no_threading);
  free(commonthreaddata);
  free(no_threading);

  method_free( &l );
  method_finalize( &l );
  
  MPI_Finalize();
  
  return 0;
}
