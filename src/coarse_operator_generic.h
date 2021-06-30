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

#ifndef COARSE_OPERATOR_PRECISION_HEADER
  #define COARSE_OPERATOR_PRECISION_HEADER

  #include "blas_vectorized.h"

  struct Thread;
  
  void coarse_operator_PRECISION_alloc( level_struct *l );
  void coarse_operator_PRECISION_free( level_struct *l );
  void coarse_operator_PRECISION_free_vectorized( operator_PRECISION_struct *op, level_struct *l );
  void coarse_operator_PRECISION_setup( vector_PRECISION *V, level_struct *l );
  void coarse_operator_PRECISION_setup_finalize( level_struct *l, struct Thread *threading );
  void coarse_operator_PRECISION_set_couplings( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_operator_PRECISION_set_self_couplings( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_operator_PRECISION_set_neighbor_couplings( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  
  void set_coarse_self_coupling_PRECISION( vector_PRECISION buffer1, vector_PRECISION buffer2,
                                           vector_PRECISION *V, const int n, level_struct *l );
  void set_coarse_neighbor_coupling_PRECISION( vector_PRECISION buffer1, vector_PRECISION buffer2,
                                               vector_PRECISION *V, const int mu, const int n, level_struct *l );

  void coarse_self_couplings_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                        operator_PRECISION_struct *op, int start, int end, level_struct *l );
  void coarse_spinwise_self_couplings_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, 
                                                 config_PRECISION clover, int length, level_struct *l );
  
  void coarse_gamma5_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start, int end, level_struct *l );
  void coarse_tau1_gamma5_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start, int end, level_struct *l );
  void apply_coarse_operator_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                        operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void g5D_apply_coarse_operator_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op,
                                            level_struct *l, struct Thread *threading );
  void apply_coarse_operator_dagger_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                               operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_block_operator_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start,
                                        schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void coarse_aggregate_self_couplings_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, 
                                                  schwarz_PRECISION_struct *s, level_struct *l );

  void coarse_aggregate_neighbor_couplings_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l );

  void set_block_diagonal_PRECISION( vector_PRECISION spin_0_1, vector_PRECISION spin_2_3, vector_PRECISION *V, const int n, config_PRECISION block, level_struct *l );

  void coarse_aggregate_block_diagonal_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, config_PRECISION block, level_struct *l );
 
  void coarse_operator_PRECISION_test_routine( level_struct *l, struct Thread *threading );
  
  // eta += D*phi, D stored columnwise
  static inline void mv_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                   const vector_PRECISION phi, const register int n ) {
    register int i, j, k=0;

    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[j] += D[k]*phi[i];
  }

  // eta -= D*phi, D stored columnwise
  static inline void nmv_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                    const vector_PRECISION phi, const register int n ) {
    register int i, j, k=0;
    
    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[j] -= D[k]*phi[i];
  }

  // eta += D^Dagger*phi, D stored columnwise
  static inline void mvh_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                    const vector_PRECISION phi, const register int n ) {
    register int i, j, k=0;    

    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[i] += conj_PRECISION(D[k])*phi[j];
  }

  // eta -= D^Dagger*phi, D stored columnwise
  static inline void nmvh_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                     const vector_PRECISION phi, const register int n ) {
    register int i, j, k=0; 

    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[i] -= conj_PRECISION(D[k])*phi[j];
  }

  // eta = D*phi, D hermitian and stored columnwise packed
  static inline void mvp_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                    const vector_PRECISION phi, const register int n ) {
    register int i, j, k;

    eta[0] = D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] = conj_PRECISION(D[k])*phi[0];
      eta[0] += D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] += D[k]*phi[i];
        eta[i] += conj_PRECISION(D[k])*phi[j];
      }
      eta[i] += D[k]*phi[i]; k++;
    }
  }

  // eta += D*phi, D hermitian and stored columnwise packed
  static inline void pmvp_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                    const vector_PRECISION phi, const register int n ) {
    register int i, j, k;

    eta[0] += D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] += conj_PRECISION(D[k])*phi[0];
      eta[0] += D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] += D[k]*phi[i];
        eta[i] += conj_PRECISION(D[k])*phi[j];
      }
      eta[i] += D[k]*phi[i]; k++;
    }
  }

  // eta += D*phi, D hermitian and stored columnwise packed
  static inline void mmvp_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                     const vector_PRECISION phi, const register int n ) {
    register int i, j, k;

    eta[0] -= D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] -= conj_PRECISION(D[k])*phi[0];
      eta[0] -= D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] -= D[k]*phi[i];
        eta[i] -= conj_PRECISION(D[k])*phi[j];
      }
      eta[i] -= D[k]*phi[i]; k++;
    }
  }

  // eta += D*phi, D anti-hermitian and stored columnwise packed
  static inline void pamvp_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                    const vector_PRECISION phi, const register int n ) {
    register int i, j, k;

    eta[0] += D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] -= conj_PRECISION(D[k])*phi[0];
      eta[0] += D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] += D[k]*phi[i];
        eta[i] -= conj_PRECISION(D[k])*phi[j];
      }
      eta[i] += D[k]*phi[i]; k++;
    }
  }
  
  // eta -= D*phi, D anti-hermitian and stored columnwise packed
  static inline void mamvp_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                    const vector_PRECISION phi, const register int n ) {
    register int i, j, k;

    eta[0] -= D[0]*phi[0];
    for ( i=1, k=1; i<n; i++ ) {
      eta[i] += conj_PRECISION(D[k])*phi[0];
      eta[0] -= D[k]*phi[i]; k++;
      for ( j=1; j<i; j++, k++ ) {
        eta[j] -= D[k]*phi[i];
        eta[i] += conj_PRECISION(D[k])*phi[j];
      }
      eta[i] -= D[k]*phi[i]; k++;
    }
  }

  static inline void coarse_self_couplings_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                                             config_PRECISION clover, int length, level_struct *l ) {
    
    int site_var = l->num_lattice_site_var,
      num_eig_vect = l->num_parent_eig_vect,
      clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2,
      clover_step_size2 = SQUARE(num_eig_vect);
    config_PRECISION clover_pt = clover;
    vector_PRECISION phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
    // U(x) = [ A B      , A=A*, D=D*, C = -B*
    //          C D ]
    // storage order: upper triangle of A, upper triangle of D, B, columnwise
    // diagonal coupling
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      while ( phi_pt < phi_end_pt ) {
        // A
        mvp_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect;//1
        phi_pt += num_eig_vect;//1
        mvp_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        // D
        eta_pt += num_eig_vect;//2
        phi_pt += num_eig_vect;//2
        clover_pt += clover_step_size1; 
        mvp_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect;//3
        phi_pt += num_eig_vect;//3
        mvp_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        // C = -B*
        eta_pt -= num_eig_vect;//2
        phi_pt -= 3*num_eig_vect;//0
        clover_pt += clover_step_size1;
        nmvh_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect;//3
        phi_pt += num_eig_vect;//1
        nmvh_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        // B
        eta_pt -= 3*num_eig_vect;//0
        phi_pt += num_eig_vect;//2
        mv_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect;//1
        phi_pt += num_eig_vect;//3
        mv_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        eta_pt += 3*num_eig_vect;//4
        phi_pt += num_eig_vect;//4
        clover_pt += clover_step_size2;
      }
    } else {
#endif


#ifdef MUMPS_ADDS

  // TODO #1 : 
  // 		sparse vals and indices : l->p_PRECISION.mumps_vals , l->p_PRECISION.mumps_Is , l->p_PRECISION.mumps_Js


    

  //To-Dos: 
// remove index k and use clover_pt ++ 


      int show_vars = 0;

//  complex_PRECISION *vals = malloc(SQUARE(site_var) * sizeof(complex_PRECISION));

      int nr_nodes = l->num_inner_lattice_sites;
      int i, j, k; // k = index in matrix
      int skip = 0;	//skip number of elements in Blockrow in large matrix (for self coupl. only skip = 0, else 8 * SQUARE(site_var))
      k = 0;
      for (j = 0; j < nr_nodes; j++){
        for (i = 0; i < SQUARE(site_var); i++, k++){
          *(l->p_PRECISION.mumps_Js +k) = j * site_var + (int)(i/site_var);	// col indices
          *(l->p_PRECISION.mumps_Is +k) = j * site_var + (i % site_var); 	// row indices
        }
        k += skip;
      }

//      printf("nr_nodes: %d, sitevar^2: %d\n", nr_nodes, SQUARE(site_var));
/*      for (k = 0; k < nr_nodes * SQUARE(site_var); k++){
        printf("I: %d, J: %d\n", *(l->p_PRECISION.mumps_Is + k), *(l->p_PRECISION.mumps_Js + k));
      }

      printf("len of Is, Js: %ld\n", nr_nodes * SQUARE(site_var));
      for (k = 0; k < nr_nodes * SQUARE(site_var); k++){
        if (*(l->p_PRECISION.mumps_Is + k) == 0) printf("found 0 in Is at k: %d\n", k);
        if (*(l->p_PRECISION.mumps_Js + k) == 0) printf("found 0 in Js at k: %d\n", k);
      }
      exit(0);
*/
      int c, r;	// col no., row no.
      for (j = 0; j < nr_nodes; j++){

	// A store column-wise
        for (k = 0, r = 0; r < num_eig_vect; r++, k++){
          for (c = 0; c < r; c++, k++){
            l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + c * site_var + r] = *(clover_pt + k);
            l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + r * site_var + c] = conj_PRECISION(*(clover_pt + k));
          }
          l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + r * site_var + r] = *(clover_pt + k);
        }


  //remove this line as well, as soon k++ is removed
        clover_pt += clover_step_size1;

	// D store column-wise
        for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
          for (c = num_eig_vect; c < r; c++, k++){
            l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + c * site_var + r] = *(clover_pt + k);
            l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + r * site_var + c] = conj_PRECISION(*(clover_pt + k));
          }
          l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + r * site_var + r] = *(clover_pt + k);
        }

  //remove this line as well, as soon k++ is removed
        clover_pt += clover_step_size1;

	// C store column-wise
        for (r = num_eig_vect, k = 0; r < 2*num_eig_vect; r++){
          for (c = 0; c < num_eig_vect; c++, k++){
            l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + (r * site_var) + c] = -1.0*(conj_PRECISION(*(clover_pt + k)));
//      l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + (r +  num_eig_vect) * site_var + c] =  -1.0*(conj_PRECISION(*(clover_pt + k)));
          }
        }

  //no clover_pt correction / change this once k++ is removed


	// B store column-wise / transposed from former storage
        for (r = 0, k = 0; r < num_eig_vect; r++){
          for (c = 0; c < num_eig_vect; c++, k++){
            l->p_PRECISION.mumps_vals[j * SQUARE(site_var) + (c * site_var) + r + num_eig_vect] = *(clover_pt + k);
          }
        }

	clover_pt += clover_step_size2;

	// here self coupl. is col-wise in mumps_vals[]

	// add skipping num for hopping terms in mumps_vals
      }  // end for loop over the blocks
      
/*

      printf("len of mumps_vals: %ld\n", nr_nodes * SQUARE(site_var));
      for (k = 0; k < nr_nodes * SQUARE(site_var); k++){
        if (*(l->p_PRECISION.mumps_vals + k) == 0) printf("found 0 in vals at k: %d\n", k);
      }
      exit(0);
*/

/*
//	old dense blas call

      complex_PRECISION my_eta[site_var];
      char N = 'T';
      complex_PRECISION one = 1, zero = 0;
      int one_i = 1;
//  extern void zgemv_(char *transA, int *m, int *n, double complex *alpha, double complex *A, int *lda, double complex *X, int *incx, double complex *beta, double  complex *Y, int *incy);
      gemv_PRECISION(&N, &site_var, &site_var, &one, &l->p_PRECISION.mumps_vals[0], &site_var, &phi[0], &one_i, &zero, &my_eta[0], &one_i);
 	 //		  T, m        , n        , alpha, A      , lda      , X         , incx  , beta , Y         , incy



*/

      clover_pt = clover;	//reset to function input
      vector_PRECISION eta_0 = eta;

#endif


      while ( phi_pt < phi_end_pt ) {


      // A
        mvp_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        clover_pt += clover_step_size1; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      // D
        mvp_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        clover_pt += clover_step_size1; phi_pt -= num_eig_vect;
      // C = -B*
        nmvh_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        phi_pt += num_eig_vect; eta_pt -= num_eig_vect;
      // B
        mv_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
        clover_pt += clover_step_size2; phi_pt += num_eig_vect; eta_pt += site_var;
      }





#ifdef MUMPS_ADDS
/* 
	//	old comparison for dense blas call
      for (i = 0; i < site_var; i++){
        printf("I:\t%d,\tDiff:\t%+f %+fi\n", i, creal(*(eta_0 +i)) - creal(*(my_eta +i)), cimag(*(eta_0 +i)) - cimag(*(my_eta +i)));
//        printf("I:\t%d,\tmy_eta:\t%+f %+fi\n", i, creal(*(my_eta +i)), cimag(*(my_eta +i)));
      }
      exit(0);
*/
#endif





    }
  }

  static inline void coarse_add_block_diagonal_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                                        config_PRECISION block, int length, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
      block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
    config_PRECISION block_pt = block;
    vector_PRECISION phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
    // U(x) = [ A 0      , A=A*, D=D* diag. excluded
    //          0 D ]
    // storage order: upper triangle of A, upper triangle of D, columnwise
    // diagonal coupling
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      while ( phi_pt < phi_end_pt ) {
        // A
        pmvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        mmvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        pmvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        mmvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      }
    } else
#endif
      while ( phi_pt < phi_end_pt ) {
        // A
        pmvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        pmvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      }
  }

  static inline void coarse_add_anti_block_diagonal_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                                               config_PRECISION block, int length, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
      block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
    config_PRECISION block_pt = block;
    vector_PRECISION phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
    // U(x) = [ A 0      , A=-A*, D=-D* diag. excluded
    //          0 D ]
    // storage order: upper triangle of A, upper triangle of D, columnwise
    // diagonal coupling
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      while ( phi_pt < phi_end_pt ) {
        // A
        pamvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        mamvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        pamvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        mamvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      }
    } else
#endif
      while ( phi_pt < phi_end_pt ) {
        // A
        pamvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
        // D
        pamvp_PRECISION( eta_pt, block_pt, phi_pt, num_eig_vect );
        block_pt += block_step_size; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
      }
  }

  static inline void coarse_add_doublet_coupling_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                                          config_PRECISION block, int length, level_struct *l ) {
    
#ifdef HAVE_TM1p1
    int num_eig_vect = l->num_parent_eig_vect,
      block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
    config_PRECISION block_pt = block;
    vector_PRECISION phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
    // U(x) = [ 0 A      , A=-A*, D=-D* diag. excluded
    //          D 0 ]
    // storage order: upper triangle of A, upper triangle of D, columnwise
    // diagonal coupling
    
    while ( phi_pt < phi_end_pt ) {
      // A
      pamvp_PRECISION( eta_pt, block_pt, phi_pt+num_eig_vect, num_eig_vect );
      pamvp_PRECISION( eta_pt+num_eig_vect, block_pt, phi_pt, num_eig_vect );
      block_pt += block_step_size; eta_pt += 2*num_eig_vect; phi_pt += 2*num_eig_vect;
      // D
      pamvp_PRECISION( eta_pt, block_pt, phi_pt+num_eig_vect, num_eig_vect );
      pamvp_PRECISION( eta_pt+num_eig_vect, block_pt, phi_pt, num_eig_vect );
      block_pt += block_step_size; eta_pt += 2*num_eig_vect; phi_pt += 2*num_eig_vect;
    }
#else
    warning0("coarse_add_doublet_coupling_PRECISION called without HAVE_TM1p1 defined.\n");
    return;
#endif
}
  
  static inline void coarse_hopp_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                            config_PRECISION D, level_struct *l ) {
  
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A  
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      // C
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//0
      D += num_eig_vect2;
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      // B
      eta -= 3*num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      // D
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//2
      D += num_eig_vect2;
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      nmv_PRECISION( eta, D, phi, num_eig_vect );
    } else {
#endif


//####################################################################
//		MY CODE
//####################################################################

      complex_PRECISION my_eta[2 * num_eig_vect];
      complex_PRECISION my_phi[2 * num_eig_vect];
      complex_PRECISION vals[4 * num_eig_vect2];
      int colInd[4 * num_eig_vect2];
      int rowPtr[4 * num_eig_vect2 +1];

      int i, j;
      for (i = 0; i < 2 * num_eig_vect; i++){
        my_eta[i] = *(eta+i);
        my_phi[i] = *(phi+i);
      }


/*
      for (i = 0; i < num_eig_vect; i++){
        printf("I:\t%d,\tDD:\t%+4.2f %+4.2fi,\town:\t%+4.2f %+4.2fi\n", i, creal(*(eta +i)), cimag(*(eta +i)), creal(my_eta[i]), cimag(my_eta[i]));
      }
      exit(0);

*/
        int k, c, r;	// index in matrix, col no., row no.
	

	// A 
        for (c = 0, k = 0; c < num_eig_vect; c++){
          for (r = 0; r < num_eig_vect; r++, k++){
            vals[r * (2*num_eig_vect) + c] = *(D + k);
          }
        }

	// C
        for (c = 0, k = 0; c < num_eig_vect; c++){
          for (r = num_eig_vect; r < 2 * num_eig_vect; r++, k++){
            vals[r * (2*num_eig_vect) + c] = *(D + num_eig_vect2 + k);
          }
        }

	// B
        for (c = num_eig_vect, k = 0; c < 2*num_eig_vect; c++){
          for (r = 0; r < num_eig_vect; r++, k++){
            vals[r * (2*num_eig_vect) + c] = *(D + 2 * num_eig_vect2 + 
k);
          }
        }

   
	// D
        for (c = num_eig_vect, k = 0; c < 2*num_eig_vect; c++){
          for (r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
            vals[r * (2*num_eig_vect) + c] = *(D + 3 * num_eig_vect2 + 
k);
          }
        }
      





      char T = 'T';
      complex_PRECISION one = 1, zero = 0, m_one = -1;
      int one_i = 1;
      int m = 2 * num_eig_vect, n = 2 * num_eig_vect;
      int lda = 2 * num_eig_vect;
//  extern void zgemv_(char *transA, int *m, int *n, double complex *alpha, double complex *A, int *lda, double complex *X, int *incx, double complex *beta, double  complex *Y, int *incy);
      gemv_PRECISION(&T, &m, &n, &m_one, &vals[0], &lda, &my_phi[0], &one_i, &one, &my_eta[0], &one_i);
  //		  T, m        , n        , alpha, A      , lda      , X         , incx  , beta , Y         , incy
   
/*
      if (show_vars == 1){
        printf("#\n");
        printf("# len of vals: \t %ld\n", sizeof(vals)/sizeof(vals[0]));	// no of els in vals
        printf("# len of colInd: \t %ld\n", sizeof(colInd)/sizeof(colInd[0]));// no of els in colInd
        printf("# len of rowPtr: \t %ld\n", sizeof(rowPtr)/sizeof(rowPtr[0]));// no of els in rowPtr
        printf("####################################################\n");
      }

      if (show_vars == 1){
        j = 0;
        for (i = 0; i < sizeof(vals)/sizeof(vals[0]); i++){
          if (vals[i] == 0) j++;
        }
        printf("found %d elements equal to 0\n", j);
      }
      */

      complex_PRECISION *eta_0 = eta;



//###########################################################
// 		ORIGINAL DDalphaAMG
//###########################################################

      // A    
      // eta -= D*phi, D stored columnwise
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      // C
      eta += num_eig_vect;
      D += num_eig_vect2;
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      // B
      phi += num_eig_vect;
      eta -= num_eig_vect;
      D += num_eig_vect2;
      nmv_PRECISION( eta, D, phi, num_eig_vect );
      // D
      eta += num_eig_vect;
      D += num_eig_vect2;
      nmv_PRECISION( eta, D, phi, num_eig_vect );


//###########################################################
// 		END OF ORIGINAL DDalphaAMG
//###########################################################

/*
      for (i = 0; i < 2 * num_eig_vect; i++){
        printf("I:\t%d,\tDiff:\t%+4.2f %+4.2fi \n", i, creal(*(eta_0 
+i)) - creal(my_eta[i]), cimag(*(eta_0 +i)) - cimag(my_eta[i]));
      }
      exit(0);
    */






#ifdef HAVE_TM1p1
    }
#endif
  }


  static inline void coarse_daggered_hopp_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                                     config_PRECISION D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A* 
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      // -C*
      eta -= num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      // -B*
      eta += num_eig_vect;//2
      phi -= 3*num_eig_vect;//0
      D += num_eig_vect2;
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      // D*
      eta -= num_eig_vect;//2
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
    } else {
#endif
      // A* 
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      // -C*
      phi += num_eig_vect;
      D += num_eig_vect2;
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      // -B*
      eta += num_eig_vect;
      phi -= num_eig_vect;
      D += num_eig_vect2;
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      // D*
      phi += num_eig_vect;
      D += num_eig_vect2;
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
#ifdef HAVE_TM1p1
    }
#endif
  }
  
  static inline void coarse_n_hopp_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                              config_PRECISION D, level_struct *l ) {
  
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A  
      mv_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      mv_PRECISION( eta, D, phi, num_eig_vect );
      // C
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//0
      D += num_eig_vect2;
      mv_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      mv_PRECISION( eta, D, phi, num_eig_vect );
      // B
      eta -= 3*num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      mv_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      mv_PRECISION( eta, D, phi, num_eig_vect );
      // D
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//2
      D += num_eig_vect2;
      mv_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      mv_PRECISION( eta, D, phi, num_eig_vect );
    } else {
#endif
      // A  
      mv_PRECISION( eta, D, phi, num_eig_vect );
      // C
      eta += num_eig_vect;
      D += num_eig_vect2;
      mv_PRECISION( eta, D, phi, num_eig_vect );
      // B
      phi += num_eig_vect;
      eta -= num_eig_vect;
      D += num_eig_vect2;
      mv_PRECISION( eta, D, phi, num_eig_vect );
      // D
      eta += num_eig_vect;
      D += num_eig_vect2;
      mv_PRECISION( eta, D, phi, num_eig_vect );
#ifdef HAVE_TM1p1
    }
#endif
  }

  static inline void coarse_n_daggered_hopp_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                                       config_PRECISION D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A* 
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      // -C*
      eta -= num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      // -B*
      eta += num_eig_vect;//2
      phi -= 3*num_eig_vect;//0
      D += num_eig_vect2;
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      // D*
      eta -= num_eig_vect;//2
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      mvh_PRECISION( eta, D, phi, num_eig_vect );
    } else {
#endif
      // A* 
      mvh_PRECISION( eta, D, phi, num_eig_vect );
      // -C*
      phi += num_eig_vect;
      D += num_eig_vect2;
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      // -B*
      eta += num_eig_vect;
      phi -= num_eig_vect;
      D += num_eig_vect2;
      nmvh_PRECISION( eta, D, phi, num_eig_vect );
      // D*
      phi += num_eig_vect;
      D += num_eig_vect2;
      mvh_PRECISION( eta, D, phi, num_eig_vect );
#ifdef HAVE_TM1p1
    }
#endif
  }

  static inline void coarse_spinwise_hopp_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, 
                                                     vector_PRECISION phi, config_PRECISION D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    // A  
    mv_PRECISION( eta1, D, phi, num_eig_vect );
    // C
    eta1 += num_eig_vect;
    D += num_eig_vect2;
    mv_PRECISION( eta1, D, phi, num_eig_vect );
    // B
    phi += num_eig_vect;
    D += num_eig_vect2;
    mv_PRECISION( eta2, D, phi, num_eig_vect );
    // D
    eta2 += num_eig_vect;
    D += num_eig_vect2;
    mv_PRECISION( eta2, D, phi, num_eig_vect );
  }


  static inline void coarse_spinwise_daggered_hopp_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2,
                                                              vector_PRECISION phi, config_PRECISION D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);  
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    // A* 
    mvh_PRECISION( eta1, D, phi, num_eig_vect );
    // -C*
    phi += num_eig_vect;
    D += num_eig_vect2;
    nmvh_PRECISION( eta2, D, phi, num_eig_vect );
    // -B*
    eta1 += num_eig_vect;
    phi -= num_eig_vect;
    D += num_eig_vect2;
    nmvh_PRECISION( eta1, D, phi, num_eig_vect );
    // D*
    eta2 += num_eig_vect;
    phi += num_eig_vect;
    D += num_eig_vect2;
    mvh_PRECISION( eta2, D, phi, num_eig_vect );
  }

  static inline void coarse_spinwise_n_hopp_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2,
                                                       vector_PRECISION phi, config_PRECISION D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    // A  
    nmv_PRECISION( eta1, D, phi, num_eig_vect );
    // C
    eta1 += num_eig_vect;
    D += num_eig_vect2;
    nmv_PRECISION( eta1, D, phi, num_eig_vect );
    // B
    phi += num_eig_vect;
    D += num_eig_vect2;
    nmv_PRECISION( eta2, D, phi, num_eig_vect );
    // D
    eta2 += num_eig_vect;
    D += num_eig_vect2;
    nmv_PRECISION( eta2, D, phi, num_eig_vect );
  }


  static inline void coarse_spinwise_n_daggered_hopp_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2,
                                                                vector_PRECISION phi, config_PRECISION D, level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);  
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    // A* 
    nmvh_PRECISION( eta1, D, phi, num_eig_vect );
    // -C*
    phi += num_eig_vect;
    D += num_eig_vect2;
    mvh_PRECISION( eta2, D, phi, num_eig_vect );
    // -B*
    eta1 += num_eig_vect;
    phi -= num_eig_vect;
    D += num_eig_vect2;
    mvh_PRECISION( eta1, D, phi, num_eig_vect );
    // D*
    eta2 += num_eig_vect;
    phi += num_eig_vect;
    D += num_eig_vect2;
    nmvh_PRECISION( eta2, D, phi, num_eig_vect );
  }

#endif
