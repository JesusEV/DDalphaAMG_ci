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
// this function computes the solution for a given LU decomposition of A stored
// in input A for a given RHS B
void coarse_scalap_solve_PRECISION(vector_PRECISION phi, vector_PRECISION Dphi,
                                   vector_PRECISION eta, int res,
                                   level_struct *l, struct Thread *threading);

// this function computes the inverse for a given Matrix A of size N x N and
// overwrites the input array with the computed LU-Decomposition
void coarse_scalap_factorize_PRECISION(level_struct *l, vector_PRECISION A,
                                       int *descA, int *ipiv,
                                       struct Thread *threading);

// this function will generate a "dense" memory layout for solving with
// scalapack
void coarse_scalap_setup_PRECISION(level_struct *l, struct Thread *threading);

void coarse_scalap_init_PRECISION(level_struct *l, struct Thread *threading);

#endif

#endif
