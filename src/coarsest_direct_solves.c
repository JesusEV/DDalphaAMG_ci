#include "main.h"


#if defined(MUMPS_ADDS) || defined(COARSE_SCALAP)

void direct_solver_init( ds_wrapper_float_struct *ds ){
    ds->mumps_vals = NULL;
    ds->mumps_Is = NULL;
    ds->mumps_Js = NULL;
    ds->mumps_rhs_loc = NULL;
    ds->mumps_irhs_loc = NULL;
    ds->mumps_SOL = NULL;
#ifdef COARSE_SCALAP
    ds->dense_vals = NULL;
    ds->dense_vals2d = NULL;
    ds->rhs2d = NULL;
    ds->desc_dense_vals = NULL;
    ds->desc_rhs = NULL;
    ds->desc_dense_vals2d = NULL;
    ds->desc_rhs2d = NULL;
    ds->ipiv = NULL;
    ds->blacs_ctxt1d = 0;
    ds->blacs_ctxt2d = 0;
    ds->myrow = 0;
    ds->myrow2d = 0;
#endif
}


void direct_solver_free( ds_wrapper_float_struct *ds, level_struct *l ){
}

#endif


