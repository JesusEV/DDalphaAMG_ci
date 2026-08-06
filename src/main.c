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

/*	
    printf("sizeof(int)      = %zu  (__SIZEOF_INT__   = %d)\n", sizeof(int),   __SIZEOF_INT__);
    printf("sizeof(long)     = %zu  (__SIZEOF_LONG__  = %d)\n", sizeof(long),  __SIZEOF_LONG__);
    printf("sizeof(void*)    = %zu\n", sizeof(void*));
    printf("sizeof(ptrdiff_t)= %zu\n", sizeof(ptrdiff_t));
    printf("__INT_MAX__      = %ld\n", (long)__INT_MAX__);
    printf("__LONG_MAX__     = %ld\n", (long)__LONG_MAX__);
    // Prüfe, ob der Compiler ein ILP64‑Makro setzt 
#ifdef __ILP64__
    puts("__ILP64__ is defined");
#else
    puts("__ILP64__ is NOT defined");
#endif
#if defined(MKL_ILP64) || defined(ILP64)
    printf("ILP64 mode is active\n");
    printf("Index type size = %zu\n", sizeof(MKL_INT));  // oder sizeof(MKL_INT)
#endif
    #ifdef ILP64
        puts("ILP64 is defined");
    #else
        puts("ILP64 is NOT defined");
    #endif

    #ifdef MKL_ILP64
        puts("MKL_ILP64 is defined");
    #else
        puts("MKL_ILP64 is NOT defined");
    #endif
*/

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

    // set up initial MG hierarchy
    method_setup( NULL, &l, &threading );

    //if ( g.my_rank == 0 ) printf("*********** p->gcrodr_PRECISION.k = %d\n", l.next_level->next_level->p_float.gcrodr_float.k);

    set_some_coarsest_level_improvs_params_for_setup( &l, &threading );

    //if ( g.my_rank == 0 ) printf("*********** p->gcrodr_PRECISION.k = %d\n", l.next_level->next_level->p_float.gcrodr_float.k);

#if defined(MUMPS_ADDS) || defined(COARSE_SCALAP)
    {
      level_struct *lx = &l;
      for (int i = 1; i<g.num_levels; i++) {
        lx = lx->next_level;
      }
      double t0 = 0,t1 = 0;
      struct Thread* threadx = &threading;
      if (!lx->idle){

	SYNC_CORES(threadx)
        START_MASTER(threadx)
	END_MASTER(threadx)
        mumps_setup_float(lx, threadx);        //setup vals, Is, Js
        START_MASTER(threadx)
        t0 = MPI_Wtime();

#ifdef MUMPS_ADDS
        //g.mumps_id.job = 4; //analyze and factorize
        g.mumps_id.job = 1; //analyze
        cmumps_c(&(g.mumps_id));
	
	if (g.on_solve) {
	        g.mumps_id.job = 2; //factorize
        	cmumps_c(&(g.mumps_id));//only factorize when on solve
	}
       
	MPI_Barrier(lx->gs_float.level_comm);
        printf0("mumps analyze + factorize done in main.c\n");
#else
	//compute LU of matrix using scalapack
      }
      if (g.on_solve){    //will be FALSE, no factorize during setup!
	    coarse_scalap_factorize_float( lx, g.ds->dense_vals2d,
		    g.ds->desc_dense_vals2d,
		    g.ds->ipiv, threadx);//only factorize when on solve
	    //printf0("scalapack factorize done in main.c\n");
      }
      if (!lx->idle){
#endif

        t1 = MPI_Wtime();
	g.coarsest_fact_time += t1 - t0;
#ifdef MUMPS_ADDS
	printf0("MUMPS analyze ");
	if (g.on_solve) printf0("and factorize ");
	printf0("time (seconds) : %f \t in main.c\n",t1-t0);
#else
	if (g.on_solve) printf0("Invert using scalapack time (seconds) : %f \t in main.c\n",t1-t0);
#endif
        END_MASTER(threadx)
        SYNC_CORES(threadx)
    
      }
    }
#endif
    MPI_Barrier(MPI_COMM_WORLD);
    printf0("starting iterative Phase\n");

    // iterative phase
    method_update( l.setup_iter, &l, &threading );

    set_some_coarsest_level_improvs_params_for_solve( &l, &threading );

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
