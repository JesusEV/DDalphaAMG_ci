#ifndef CDS_PRECISION_HEADER
  #define CDS_PRECISION_HEADER

  #define ICNTL(I) icntl[(I) -1]	//macro according to docu //bridges from fortran indices to c

// this function will set up the data format for mumps / csr
void mumps_setup_PRECISION(level_struct *l, struct Thread *threading);

// this function will set / reset some important variables. It is called in
// top_level.c, once the setup phase is done and right before the solve phase
void direct_solves_set_reset_PRECISION(level_struct *l,
                                       struct Thread *threading);

#ifdef MUMPS_ADDS
// this function will do all the necessary handling of data for the solve call.
// e.g. distributing the calculated solution to all processes
void mumps_solve_PRECISION(vector_PRECISION phi, vector_PRECISION Dphi,
                           vector_PRECISION eta, int res, level_struct *l,
                           struct Thread *threading);

// this function will initialize the cmumps instance, set control parameter,
// general values and link arrays
void mumps_init_PRECISION(gmres_PRECISION_struct *p, int mumps_n, int nnz_loc,
                          int rhs_len, level_struct *l, Thread *threading);
#endif

#ifdef COARSE_SCALAP
// Direct solver on the odd-even Schur complement S = D_ee - D_eo D_oo^{-1} D_oe

// solves D x = b on the coarsest level: b_e <- b_e - D_eo D_oo^{-1} b_o,
// x_e = S^{-1} b_e (dense, scalapack), x_o = D_oo^{-1} (b_o - D_oe x_e).
// b is overwritten. Collective over all processes.
void coarse_scalap_solve_oddeven_PRECISION(vector_PRECISION x, vector_PRECISION b,
                                           level_struct *l, struct Thread *threading);

// x_e = S^{-1} b_e on the even parts of x and b. Collective over all processes.
void coarse_scalap_apply_Sinv_PRECISION(vector_PRECISION x, vector_PRECISION b,
                                        level_struct *l, struct Thread *threading);

// computes S = D_ee - D_eo (D_oo^{-1} D_oe) with pgemm and overwrites it with
// its explicit inverse (pgetrf + pgetri). Collective over all processes.
void coarse_scalap_factorize_PRECISION(level_struct *l, struct Thread *threading);

// sorts the locally assembled rows of D into the 1D blocks D_ee, D_eo, D_oe,
// D_zz and computes D_oo^{-1} D_oe locally (no communication)
void coarse_scalap_setup_PRECISION(level_struct *l, struct Thread *threading);

// BLACS grids, descriptors, one-time allocation of the 2D arrays and workspace
void coarse_scalap_init_PRECISION(level_struct *l, struct Thread *threading);

// one-time test: ||b - D x||/||b|| and ||x_e - S^{-1} S x_e||/||x_e||
void coarse_scalap_selftest_PRECISION(level_struct *l);

// moves D_ee, D_eo and D_oo^{-1} D_oe from the 1D to the 2D block cyclic pattern
void scalap_1d_2d_A_PRECISION(level_struct *l, struct Thread *threading);

// moves the even part of vec from the 1D pattern to rhs2d (2D block cyclic)
void scalap_1d_2d_vec_PRECISION(vector_PRECISION vec, level_struct *l,
                                struct Thread *threading);

// moves sol2d (2D block cyclic) to the even part of vec (1D pattern)
void scalap_2d_1d_vec_PRECISION(vector_PRECISION vec, level_struct *l,
                                struct Thread *threading);

#endif

#endif
