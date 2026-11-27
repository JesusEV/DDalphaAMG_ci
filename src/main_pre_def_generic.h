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

#ifndef MAIN_PRE_DEF_PRECISION_HEADER
  #define MAIN_PRE_DEF_PRECISION_HEADER
  
  typedef PRECISION complex complex_PRECISION;
  typedef PRECISION complex *config_PRECISION;
  typedef PRECISION complex *vector_PRECISION;

  typedef struct {
    int length[8], *boundary_table[8], max_length[4],
        comm_start[8], in_use[8], offset, comm,
        num_even_boundary_sites[8], num_odd_boundary_sites[8],
        num_boundary_sites[8];
    vector_PRECISION buffer[8];
    MPI_Request sreqs[8], rreqs[8];
  } comm_PRECISION_struct;
  
  typedef struct {
    int ilde, dist_local_lattice[4], dist_inner_lattice_sites,
        *permutation, *gather_list, gather_list_length;
    vector_PRECISION buffer, transfer_buffer;
    MPI_Request *reqs;
    MPI_Group level_comm_group;
    MPI_Comm level_comm;
  } gathering_PRECISION_struct;
  
  typedef struct {
    double m0;
    config_PRECISION D, clover, clover_oo_inv;
    config_PRECISION odd_proj; //identity on the odd sites
    int oe_offset, self_coupling, num_even_sites, num_odd_sites,
        *index_table, *neighbor_table, *translation_table, table_dim[4],
        *backward_neighbor_table,
        table_mod_dim[4], *config_boundary_table[4];
    vector_PRECISION *buffer, prnT, prnZ, prnY, prnX, prpT, prpZ, prpY, prpX;
    comm_PRECISION_struct c;
    OPERATOR_TYPE_PRECISION *D_vectorized;
    OPERATOR_TYPE_PRECISION *D_transformed_vectorized;
    OPERATOR_TYPE_PRECISION *clover_vectorized;
    OPERATOR_TYPE_PRECISION *clover_oo_inv_vectorized;
#ifdef HAVE_TM
    double mu, mu_odd_shift, mu_even_shift;
    config_PRECISION tm_term;
#endif
#ifdef HAVE_TM1p1
    double epsbar, epsbar_ig5_odd_shift, epsbar_ig5_even_shift;
    config_PRECISION epsbar_term, clover_doublet_oo_inv;
    OPERATOR_TYPE_PRECISION *clover_doublet_vectorized;
    OPERATOR_TYPE_PRECISION *clover_doublet_oo_inv_vectorized;
#endif
  } operator_PRECISION_struct;

#if defined(POLYPREC) || defined(GCRODR)
  typedef struct
  {
#ifdef COARSE_SCALAP
    lapack_int N, nrhs, lda, ldb, info;
    lapack_int *ipiv;
#else
    int N, nrhs, lda, ldb, info;
    int *ipiv;
#endif
    vector_PRECISION x, b;
    complex_PRECISION *Hcc;  

    void (*dirctslvr_PRECISION)();

  } dirctslvr_PRECISION_struct;
#endif

#if defined(GCRODR) || defined(POLYPREC)
  // this is both eigensolver and generalized eigensolver
  typedef struct {
    char jobvl, jobvr;

#ifdef COARSE_SCALAP
    lapack_int N, lda, ldb, ldvl, ldvr, info, qr_m, qr_n, qr_lda, qr_k;
#else
    int N, lda, ldb, ldvl, ldvr, info, qr_m, qr_n, qr_lda, qr_k;
#endif
    int *ordr_idxs;

    complex_PRECISION *ordr_keyscpy, *qr_tau;
    vector_PRECISION vl, vr, w, beta, A, B;

    complex_PRECISION **qr_QR, **qr_Q, **qr_R, **qr_Rinv;
    complex_PRECISION **Hc;

    void (*eigslvr_PRECISION)();
    void (*gen_eigslvr_PRECISION)();
  } eigslvr_PRECISION_struct;
#endif

#ifdef GCRODR
  typedef struct {
    int i, k, CU_usable, syst_size, finish, orth_against_Ck, update_CU, recompute_DPCk_poly, recompute_DPCk_plain, upd_ctr;

    PRECISION b_norm, norm_r0;

    vector_PRECISION *Pk, *C, *Cc, *U, *Yk, *hatZ, *hatW;
//#ifdef BLOCK_JACOBI
#if 0
    vector_PRECISION r_aux;
#endif
    // Gc is used to copy G
    complex_PRECISION *lsp_x, *lsp_diag_G, **lsp_H;
    complex_PRECISION **gev_A, **gev_B, **Bbuff, **QR, **Q, **R, **Rinv, **ort_B, **G, **Gc;

    eigslvr_PRECISION_struct eigslvr;

#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    vector_PRECISION *PC, *DPC;
#endif
  } gcrodr_PRECISION_struct;
#endif

#ifdef POLYPREC
  typedef struct
  {
    int update_lejas;
    int d_poly;
    int syst_size;
      
    complex_PRECISION **Hc;
    complex_PRECISION *Hcc;
    complex_PRECISION **L;
    complex_PRECISION *col_prods;
    vector_PRECISION h_ritz;
    vector_PRECISION lejas;
    vector_PRECISION random_rhs;
    vector_PRECISION accum_prod, product, temp, xtmp;

    void (*preconditioner)();
    void (*preconditioner_bare)();

    eigslvr_PRECISION_struct eigslvr;
    dirctslvr_PRECISION_struct dirctslvr;
  } polyprec_PRECISION_struct;
#endif

//#ifdef BLOCK_JACOBI
#if 0
  typedef struct {
    vector_PRECISION x, b, r, w, *V, *Z;
    complex_PRECISION **H, *y, *gamma, *c, *s;
    config_PRECISION *D, *clover;
    operator_PRECISION_struct *op;
    PRECISION tol;
    int num_restart, restart_length, timing, print, kind,
      initial_guess_zero, layout, v_start, v_end;
    long int total_storage;
    void (*eval_operator)();

    polyprec_PRECISION_struct polyprec_PRECISION;
  } local_gmres_PRECISION_struct;

  typedef struct {
    int BJ_usable, syst_size;
    vector_PRECISION b_backup;
    vector_PRECISION xtmp;
    local_gmres_PRECISION_struct local_p;

    // for direct solves
    OPERATOR_TYPE_PRECISION* bj_op_inv_vectorized;
    OPERATOR_TYPE_PRECISION* bj_op_vectorized;
    OPERATOR_TYPE_PRECISION* bj_doublet_op_inv_vectorized;
    OPERATOR_TYPE_PRECISION* bj_doublet_op_vectorized;

    //vector_PRECISION xxxtmp[4];

  } block_jacobi_PRECISION_struct;
#endif

  typedef struct {
    vector_PRECISION x, b, r, w, *V, *Z, richardson_w, richardson_r;
    complex_PRECISION **H, *y, *gamma, *c, *s;
    config_PRECISION *D, *clover;
    operator_PRECISION_struct *op;
    PRECISION tol, richardson_omega;
    int num_restart, restart_length, timing, print, kind,
      initial_guess_zero, layout, v_start, v_end, richardson_update_omega,
      richardson_vl;
    long int total_storage;
    void (*preconditioner)();
    void (*eval_operator)();
    
#ifdef GCRODR
    gcrodr_PRECISION_struct gcrodr_PRECISION;
    vector_PRECISION rhs_bk;
    int was_there_stagnation;
#endif
#ifdef POLYPREC
    polyprec_PRECISION_struct polyprec_PRECISION;
#endif
//#ifdef BLOCK_JACOBI
#if 0
    block_jacobi_PRECISION_struct block_jacobi_PRECISION;
#endif
#if defined(SINGLE_ALLREDUCE_ARNOLDI) && defined(PIPELINED_ARNOLDI)
    int syst_size;
    vector_PRECISION *Va, *Za;
#endif

#ifdef PERS_COMMS
    vector_PRECISION* pers_comms_ins;
    vector_PRECISION* pers_comms_outs;
#endif
  } gmres_PRECISION_struct;

#if defined(MUMPS_ADDS) || defined(COARSE_SCALAP)
  typedef struct {
    vector_PRECISION mumps_vals;
    int *mumps_Is, *mumps_Js;
    vector_PRECISION mumps_rhs_loc;
    int *mumps_irhs_loc;
    vector_PRECISION mumps_SOL;
#ifdef COARSE_SCALAP
    // Direct solver on the odd-even Schur complement
    //   S = D_ee - D_eo D_oo^{-1} D_oe
    // Since every local coarsest lattice has as many even as odd sites, all
    // four blocks (and S) have the same global size N_h x N_h with N_h = N/2,
    // so one descriptor per layout (1D / 2D) serves all of them.

    // 1D layout (DDalphaAMG native, process r holds the block rows of its own sites)
    vector_PRECISION S1d;           // initially D_ee (local even rows, all even columns)
    vector_PRECISION Deo1d;         // D_eo           (local even rows, all odd columns)
    vector_PRECISION DooinvDoe1d;   // D_oo^{-1} D_oe (local odd rows, all even columns)
    vector_PRECISION Doo_blocks;    // diagonal blocks D_zz of the local odd sites (column major)
    // 2D block cyclic layout (computations)
    vector_PRECISION S2d;           // S, after coarse_scalap_factorize: S^{-1}
    vector_PRECISION Deo2d;
    vector_PRECISION DooinvDoe2d;
    vector_PRECISION rhs2d, sol2d;  // even part of rhs / solution
    vector_PRECISION work;          // workspace for pgetri
    lapack_int *iwork;              // workspace for pgetri
    lapack_int *ipiv;
    lapack_int desc_mat1d[9], desc_vec1d[9];
    lapack_int desc_mat2d[9], desc_vec2d[9];
    lapack_int N_h, bs_h;           // global size of S and 1D block size (= rows per process)
    long len_mat1d, len_doo, len_mat2d, len_vec2d, len_ipiv, lwork, liwork; // allocation lengths
    lapack_int blacs_ctxt1d;
    lapack_int blacs_ctxt2d;
    lapack_int myrow, myrow2d, mycol2d;
#endif
  } ds_wrapper_PRECISION_struct; // wrapper for direct solver (for coarsest
                                 // level)
#endif

  typedef struct {
    operator_PRECISION_struct op;
    vector_PRECISION buf1, buf2, buf3, buf4, buf5;
    vector_PRECISION oe_buf[4];
    vector_PRECISION local_minres_buffer[3];
    int block_oe_offset, *index[4], dir_length[4], num_blocks, num_colors,
        dir_length_even[4], dir_length_odd[4], *oe_index[4],
        num_block_even_sites, num_block_odd_sites, num_aggregates,
        block_vector_size, num_block_sites, block_boundary_length[9],
        **block_list, *block_list_length;
    block_struct *block;
  } schwarz_PRECISION_struct;
  
  typedef struct {
    int num_agg, *agg_index[4], agg_length[4], *agg_boundary_index[4],
        *agg_boundary_neighbor[4], agg_boundary_length[4], num_bootstrap_vect;
    vector_PRECISION *test_vector, *interpolation, *bootstrap_vector, tmp;
    complex_PRECISION *operator, *eigenvalues, *bootstrap_eigenvalues;
  } interpolation_PRECISION_struct;
  
  typedef struct {
    double time[_NUM_PROF];
    double flop[_NUM_PROF];
    double count[_NUM_PROF];
    char name[_NUM_PROF][50];
  } profiling_PRECISION_struct;
  
  #ifdef PROFILING
    #define PROF_PRECISION_START_UNTHREADED( TYPE ) do{ l->prof_PRECISION.time[TYPE] -= MPI_Wtime(); }while(0)
    #define PROF_PRECISION_START_THREADED( TYPE, threading ) do{ if(threading->core + threading->thread == 0) l->prof_PRECISION.time[TYPE] -= MPI_Wtime(); }while(0)
  #else
    #define PROF_PRECISION_START_UNTHREADED( TYPE )
    #define PROF_PRECISION_START_THREADED( TYPE, threading )
  #endif
  
  #ifdef PROFILING
    #define PROF_PRECISION_STOP_UNTHREADED( TYPE, COUNT ) do{ l->prof_PRECISION.time[TYPE] += MPI_Wtime(); \
    l->prof_PRECISION.count[TYPE] += COUNT; }while(0)
    #define PROF_PRECISION_STOP_THREADED( TYPE, COUNT, threading ) do{ if(threading->core + threading->thread == 0) { l->prof_PRECISION.time[TYPE] += MPI_Wtime(); \
    l->prof_PRECISION.count[TYPE] += COUNT; } }while(0)
  #else
    #define PROF_PRECISION_STOP_UNTHREADED( TYPE, COUNT )
    #define PROF_PRECISION_STOP_THREADED( TYPE, COUNT, threading )
  #endif

  #define GET_MACRO2(_1,_2,NAME,...) NAME
  #define GET_MACRO3(_1,_2,_3,NAME,...) NAME

  #define PROF_PRECISION_START(...) GET_MACRO2(__VA_ARGS__, PROF_PRECISION_START_THREADED, PROF_PRECISION_START_UNTHREADED, padding)(__VA_ARGS__)
  #define PROF_PRECISION_STOP(...) GET_MACRO3(__VA_ARGS__, PROF_PRECISION_STOP_THREADED, PROF_PRECISION_STOP_UNTHREADED, padding)(__VA_ARGS__)

#endif
