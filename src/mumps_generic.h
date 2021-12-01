#ifndef MUMPS_PRECISION_HEADER
  #define MUMPS_PRECISION_HEADER

  #define ICNTL(I) icntl[(I) -1]	//macro according to docu //bridges from fortran indices to c

  void mumps_setup_PRECISION(level_struct *l, struct Thread *threading);

  void mumps_solve_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                              int res, level_struct *l, struct Thread *threading );


#endif
