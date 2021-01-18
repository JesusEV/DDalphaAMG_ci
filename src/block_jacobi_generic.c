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

#ifdef BLOCK_JACOBI


  // declarations of aux functions
  // TODO


  // -------------------------------------------------------------


  // main functions here -- IMPORTANT : all Block Jacobi functions are, for now,
  //			                exclusively per-process operations

  void block_jacobi_apply_PRECISION( vector_PRECISION out, vector_PRECISION in, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ){

    START_MASTER(threading)

    if ( p->block_jacobi_PRECISION.BJ_usable==0 ) { return; }
    else if ( p->block_jacobi_PRECISION.BJ_usable==1 ) {

      //vector_PRECISION in_x=0;

      if ( out==in ) {
        // TODO :
        // 		1. assign in_x to the BJ buffer
        //		2. copy <in> into <in_x>
      } //else { in_x = in; }

      // TODO : apply Block Jacobi, with input <in_x> and output <output>

    }

    END_MASTER(threading)

  }

  void block_jacobi_restore_from_buffer_PRECISION( vector_PRECISION out, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

    START_MASTER(threading)

    //vector_PRECISION out_x=0;

    // TODO :
    // 		1. assign out_x to the BJ buffer (the same buffer used by block_jacobi_apply_PRECISION(...))
    //		2. copy <out_x> into <out>

    END_MASTER(threading)

  }

#endif
