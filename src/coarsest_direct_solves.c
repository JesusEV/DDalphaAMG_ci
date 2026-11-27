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
    ds->S1d = NULL;
    ds->Deo1d = NULL;
    ds->DooinvDoe1d = NULL;
    ds->Doo_blocks = NULL;
    ds->S2d = NULL;
    ds->Deo2d = NULL;
    ds->DooinvDoe2d = NULL;
    ds->rhs2d = NULL;
    ds->sol2d = NULL;
    ds->work = NULL;
    ds->iwork = NULL;
    ds->ipiv = NULL;
    for (int i = 0; i < 9; i++) {
        ds->desc_mat1d[i] = 0; ds->desc_vec1d[i] = 0;
        ds->desc_mat2d[i] = 0; ds->desc_vec2d[i] = 0;
    }
    ds->N_h = 0;
    ds->bs_h = 0;
    ds->len_mat1d = 0;
    ds->len_doo = 0;
    ds->len_mat2d = 0;
    ds->len_vec2d = 0;
    ds->len_ipiv = 0;
    ds->lwork = 0;
    ds->liwork = 0;
    ds->blacs_ctxt1d = 0;
    ds->blacs_ctxt2d = 0;
    ds->myrow = -1;
    ds->myrow2d = -1;
    ds->mycol2d = -1;
#endif
}


void direct_solver_free( ds_wrapper_float_struct *ds, level_struct *l ){
}

#endif

