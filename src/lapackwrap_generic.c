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

#if defined(GCRODR) || defined(POLYPREC)


// inversion of a triangular matrix
void inv_tri_PRECISION(eigslvr_PRECISION_struct* eigen_struct)
{

  //#define trtri_double LAPACKE_ztrtri
  trtri_PRECISION( LAPACK_COL_MAJOR, 'U', 'N', eigen_struct->qr_k, eigen_struct->qr_R[0],
                   eigen_struct->qr_k);
}

void qr_PRECISION(eigslvr_PRECISION_struct* eigen_struct)
{

  //#define geqr2_double LAPACKE_zgeqr2
  geqr2_PRECISION( LAPACK_COL_MAJOR, eigen_struct->qr_m, eigen_struct->qr_n, eigen_struct->qr_QR[0],
                   eigen_struct->qr_lda, eigen_struct->qr_tau );

}

void q_from_qr_PRECISION(eigslvr_PRECISION_struct* eigen_struct)
{

  //#define ungqr_double LAPACKE_zungqr
  ungqr_PRECISION( LAPACK_COL_MAJOR, eigen_struct->qr_m, eigen_struct->qr_n, eigen_struct->qr_k,
                   eigen_struct->qr_QR[0], eigen_struct->qr_lda, eigen_struct->qr_tau );

}

void gen_eigslvr_PRECISION(eigslvr_PRECISION_struct* eigen_struct)
{

  eigen_struct->info = ggev_PRECISION( LAPACK_COL_MAJOR, eigen_struct->jobvl, eigen_struct->jobvr,
                                       eigen_struct->N, eigen_struct->A, eigen_struct->lda,
                                       eigen_struct->B, eigen_struct->ldb, eigen_struct->w, eigen_struct->beta,
                                       eigen_struct->vl, eigen_struct->ldvl, eigen_struct->vr, eigen_struct->ldvr );
}

void eigslvr_PRECISION(eigslvr_PRECISION_struct* eigen_struct)
{

  eigen_struct->info = geev_PRECISION( LAPACK_ROW_MAJOR, eigen_struct->jobvl, eigen_struct->jobvr,
                                      eigen_struct->N, eigen_struct->A, eigen_struct->lda, eigen_struct->w,
                                      eigen_struct->vl, eigen_struct->ldvl, eigen_struct->vr, eigen_struct->ldvr );       
}

void dirctslvr_PRECISION(dirctslvr_PRECISION_struct* dirctslvr)
{

    memcpy( dirctslvr->x, dirctslvr->b, sizeof(complex_PRECISION)*(dirctslvr->N) );  

    dirctslvr->info = gesv_PRECISION( LAPACK_COL_MAJOR, dirctslvr->N, dirctslvr->nrhs, dirctslvr->Hcc, 
                                    dirctslvr->lda, dirctslvr->ipiv, dirctslvr->x, dirctslvr->ldb ); 
}

#endif
