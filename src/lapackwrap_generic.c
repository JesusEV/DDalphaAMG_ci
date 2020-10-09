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

#endif
