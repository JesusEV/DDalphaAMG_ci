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

#if defined(MUMPS_ADDS) || defined(COARSE_SCALAP)

#include "coarsest_direct_solves_PRECISION.h"

#ifdef COARSE_SCALAP

void blacs_get_(const lapack_int*, const lapack_int*, lapack_int*);

void blacs_pinfo_(lapack_int*, lapack_int*);
void blacs_gridinit_(lapack_int*, const char*, const lapack_int*, const lapack_int*);
void blacs_gridmap_(lapack_int *ictxt, lapack_int *usermap, lapack_int *ldumap, lapack_int *nprow, lapack_int *npcol);
void blacs_gridinfo_(const lapack_int*, lapack_int*, lapack_int*, lapack_int*, lapack_int*);
void descinit_(lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, lapack_int*);
lapack_int numroc_(const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*);

void pgesv_PRECISION(const lapack_int*, const lapack_int*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, lapack_int*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, lapack_int* );
void pgetrf_PRECISION(const lapack_int*, const lapack_int*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, lapack_int*, lapack_int* );
//		void(const lapack_int *, const lapack_int *, _Complex float *, const lapack_int *, const lapack_int *, const lapack_int *, lapack_int *, lapack_int *)
void pgetrs_PRECISION(const char*, const lapack_int*, const lapack_int*, const complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, lapack_int* );
//		 void(const char *, const int *, const int *, const _Complex float *, const int *, const int *, const int *, const int *, _Complex float *, const int *, const int *, const int *, int *)
void pgetri_PRECISION( const lapack_int*, complex_PRECISION*, const lapack_int*, const lapack_int*, lapack_int*, lapack_int*, complex PRECISION*, const lapack_int*, const lapack_int*, lapack_int*, lapack_int* );
// alpha and beta are complex for pcgemv_/pzgemv_
void pgemv_PRECISION(char*, lapack_int*, lapack_int*, complex_PRECISION*, complex_PRECISION*, lapack_int*, lapack_int*, lapack_int*, complex_PRECISION*, lapack_int*, lapack_int*, lapack_int*, lapack_int*, complex_PRECISION*, complex_PRECISION*, lapack_int*, lapack_int*, lapack_int*, lapack_int*);
void pgemm_PRECISION(char*, char*, lapack_int*, lapack_int*, lapack_int*, complex_PRECISION*, complex_PRECISION*, lapack_int*, lapack_int*, lapack_int*, complex_PRECISION*, lapack_int*, lapack_int*, lapack_int*, complex_PRECISION*, complex_PRECISION*, lapack_int*, lapack_int*, lapack_int*);
void pgemr2d_PRECISION( const lapack_int*, const lapack_int*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, lapack_int*);

void paxpy_PRECISION( const lapack_int*, complex_PRECISION*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*);
void pnrm2_PRECISION( const lapack_int*, PRECISION*, complex_PRECISION*, const lapack_int*, const lapack_int*, const lapack_int*, const lapack_int*);
#endif

void mumps_setup_PRECISION(level_struct *l, struct Thread *threading){
  /* Setting Mumps data-storage
   *
   * l will be coarsest level
   */

  SYNC_CORES(threading)
  START_MASTER(threading)
  //printf0("starting setup!\n");
  
  int lrank = 0; //local rank in level-comm
  MPI_Comm_rank(l->gs_PRECISION.level_comm, &lrank);

  // variables for timing
  double t0,t1;
  t0 = MPI_Wtime();
  
  // use variable to have a better overview
  gmres_PRECISION_struct* px = &(l->p_PRECISION);
  operator_PRECISION_struct* op = px->op;
  config_PRECISION clover_pt = l->p_PRECISION.op->clover;

  int num_eig_vect = l->num_parent_eig_vect;
  int site_var = l->num_lattice_site_var,    // contains number of vector elements per lattice site
      clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2,	//contains number of elements of
          // clover part, which is stored in triangular form / diagonal clover part
      clover_step_size2 = SQUARE(num_eig_vect);  //contains number of elements of clover part, which
          //  ist stored as dense matrix / "off-diagonal" clover part

  int nr_nodes = l->num_inner_lattice_sites;
  int i, j, k; // k = index in matrix
  int c, r;    // col no., row no.
  int skip = 8 * SQUARE(site_var); //skip number of elements in Blockrow in large matrix (for self
      //  coupl. only skip = 0, else 8 * SQUARE(site_var)), (8 = T+, T-, X+, X-,...)


  //TODO use correct rank. eg. rank on coarsest level
  // indices/positions to copy from/to
  int j_start = nr_nodes * lrank * site_var, i_start = nr_nodes * lrank * site_var; 
  
  // putting indices/positions in matrix for SELF COUPLING = CLOVER
  for (j = 0, k = 0; j < nr_nodes; j++){
    for (i = 0; i < SQUARE(site_var); i++, k++){
      *(g.ds->mumps_Is +k) = i_start + j * site_var + (int)(i/site_var);	// col indices
      *(g.ds->mumps_Js +k) = j_start + j * site_var + (i % site_var); 	// row indices
    }
    k += skip;
  }
  //printf0("inital counting done!\n");

  // A B
  // C D    all A, B, C, D stored columnwise
  // putting values to mumps formatting for SELF COUPLING = CLOVER
  for (j = 0; j < nr_nodes; j++){
    // A
    for (k = 0, r = 0; r < num_eig_vect; r++, k++){
      for (c = 0; c < r; c++, k++){
        g.ds->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] = *(clover_pt + k); //clover is
            // triangular store in DDalphaAMG
        g.ds->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] = conj_PRECISION(*(clover_pt + k));
      }
      g.ds->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] = *(clover_pt + k); // diagonal element
    }
    clover_pt += clover_step_size1; // bend pointer to next piece of memory/clover part
    // D
    for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
      for (c = num_eig_vect; c < r; c++, k++){
        g.ds->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] = *(clover_pt + k);
        g.ds->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] = conj_PRECISION(*(clover_pt + k));
      }
      g.ds->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] = *(clover_pt + k);
    }
    clover_pt += clover_step_size1;
    // C
    for (r = num_eig_vect, k = 0; r < 2*num_eig_vect; r++){
      for (c = 0; c < num_eig_vect; c++, k++){
        g.ds->mumps_vals[j*9*SQUARE(site_var) + (r * site_var) + c] = -1.0*(conj_PRECISION(*(clover_pt + k)));
      }
    }

    // no clover_pt correction, clover is hermitian, use same data as for C 
    // B store column-wise / transposed from former storage
    for (r = 0, k = 0; r < num_eig_vect; r++){
      for (c = 0; c < num_eig_vect; c++, k++){
        g.ds->mumps_vals[j*9*SQUARE(site_var) + (c * site_var) + r + num_eig_vect] = *(clover_pt + k);
      }
    }
    clover_pt += clover_step_size2; // bend pointer to next clover for next lattice site
  }

  //printf0("clover part done!\n");

#ifdef HAVE_TM
  // twisted mass-term:
  // correction of A 0
  //               0 D 
  int block_step_size = (num_eig_vect * (num_eig_vect+1))/2;
  config_PRECISION tm_block_pt = op->tm_term;  // pointer to tm_term beginning

  for (j = 0; j < nr_nodes; j++){
    // A
    for (k = 0, r = 0; r < num_eig_vect; r++, k++){
      for (c = 0; c < r; c++, k++){
        g.ds->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] += *(tm_block_pt + k);
        g.ds->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] -= conj_PRECISION(*(tm_block_pt + k));
      }
      g.ds->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] += *(tm_block_pt + k);
    }
    tm_block_pt += block_step_size;

    // D
    for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
      for (c = num_eig_vect; c < r; c++, k++){
        g.ds->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] += *(tm_block_pt + k);
        g.ds->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] -= conj_PRECISION(*(tm_block_pt + k));
      }
      g.ds->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] += *(tm_block_pt + k);
    }
    tm_block_pt += block_step_size;
  }
 // printf0("twisted mass part done!\n");
#endif

  // hopping-term
  // memory for mumps will look like: 
  /*	vals = [[self_coupling of site 1][T-_coupling site 1][T+_coupling site 1][Z-_coupling site 1][Z+_coupling site 1] .... [X-_coupling site 1][X+_coupling site 1]
  [self_coup site 2][T-_coup site 2]....[X+_coup site N]]
  each of the inner [] contain num_link_var elements -> to store 1 block in matrix (entire coupling of one site) we need 9 * num_link_var elements
  vals =   [[self, T-, T+, Z-, Z+, Y-, Y+, X-, X+][self, T-, T+, Z-, Z+, Y-, Y+, X-, X+]....]
  */

  int index, // will be used to access the neighbour_table
      num_4link_var=4*4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      num_link_var=4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      start=0;

  int core_start = start; // start and end lattice site number for each process. 
  int core_end = start+nr_nodes;
  int comm_nr[4] = {0, 0, 0, 0}; // will hold the number elements to communicate in the corresponding direction
  int dir;
  int node;

  // count the number of elements to communicate in the corresponding direction:
  for (dir = T; dir <= X; dir++){
    for (node = core_start; node < core_end; node++){
      if (op->neighbor_table[5*node+1+dir] >= nr_nodes){
        comm_nr[dir]++;
      }
    }
  }


  
  int global_comm_size;
  MPI_Comm_size(MPI_COMM_WORLD, &global_comm_size);

  int *glob_ranks = NULL;
  int *loc_ranks = NULL; //will contain global ranks and corresponding local ranks

//TODO: fix threading in here!
//  START_MASTER(threading)
  MALLOC( glob_ranks, int, global_comm_size);
  MALLOC( loc_ranks, int, global_comm_size);
//  END_MASTER(threading)
//  SYNC_CORES(threading)
 
  for (i = 0; i<global_comm_size; i++) { glob_ranks[i]= i; loc_ranks[i] = -1;}

  MPI_Group world_group;
  MPI_Group local_group;

  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Comm_group(l->gs_PRECISION.level_comm, &local_group);

  MPI_Group_translate_ranks(world_group, global_comm_size, glob_ranks, local_group, loc_ranks);
  //translates from global to local ranks

  //allocate memory for buffers:
  int *buff_i_send[4], *buff_i_recv[4]; // will contain node number in receiver processors domain
  complex_PRECISION *buff_d_send[4], *buff_d_recv[4]; // will contain mu- coupling

  int buffer_i_pt, buffer_d_pt; // hold the sweeping index in buffer
  MPI_Request req;
  MPI_Status s;

  i_start = lrank * l->num_inner_lattice_sites * site_var; 
  j_start = lrank * l->num_inner_lattice_sites * site_var; //contains own global row and col. start indices
  int neighbors_j_start; //contains column start index of neighboring process
  
  int *boundary_table; // boundary table holds the neighboring lattice sites in each direction? 
  int bt_index;
  int num_site_var=site_var;

  int li = 0; //lex index for site

  //printf0("allocs for hopping done!\n");
  
  /* FINDME */
  for (dir = T; dir <= X; dir++){
    boundary_table = op->c.boundary_table[2*dir];
    buffer_i_pt = 0;

    buff_i_send[dir] = NULL;
    buff_i_recv[dir] = NULL;
    buff_d_send[dir] = NULL;
    buff_d_recv[dir] = NULL;
 
    // allocating and initialization of buffers
    // TODO: fix threading
  //  START_MASTER(threading)
    MALLOC(buff_i_send[dir], int, 2 * comm_nr[dir]);
    MALLOC(buff_i_recv[dir], int, 2 * comm_nr[dir]);
    MALLOC(buff_d_send[dir], complex_PRECISION, num_link_var * comm_nr[dir]);
    MALLOC(buff_d_recv[dir], complex_PRECISION, num_link_var * comm_nr[dir]);
//    END_MASTER(threading)
//    SYNC_CORES(threading)
    memset(buff_i_send[dir], 0, 2 * comm_nr[dir] * sizeof(int));
    memset(buff_i_recv[dir], 0, 2 * comm_nr[dir] * sizeof(int));
    memset(buff_d_send[dir], 0, num_link_var * comm_nr[dir] * sizeof(complex_PRECISION));
    memset(buff_d_recv[dir], 0, num_link_var * comm_nr[dir] * sizeof(complex_PRECISION));


    bt_index = 0;
    // compute column index for coupling to a lattice site on a neighboring process
    // TODO: when odd-even enabled -> change this? 
    neighbors_j_start = loc_ranks[l->neighbor_rank[2*dir]] * l->num_inner_lattice_sites * site_var;

    for (node = core_start; node < core_end; node ++){
      index = 5 * node; // neighbor table will contain site numbers in chunks of 5 for each lattice site: [my_site_number, T neighbor, Z neighbor, Y neighbor, X neighbor]
    
      // make mu+ couplings as usual (Values + Row indices aka. Is)
      // A
      for (k = 0; k < SQUARE(num_site_var/2); k ++){ 
        //find correct block row       skip self coupl., find pos of mu+ coupl. ("2*mu +1" due to structure of vals[self, T-, T+, Z-, Z+...]
        *(g.ds->mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir + 1)*num_link_var + k) = 
			-1.0 * *(op->D + 	num_4link_var*op->neighbor_table[index] +
			dir*num_link_var + k); 
	//					find correct block row				start of mu- coupling
        *(g.ds->mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) = 
			i_start +	num_site_var * op->neighbor_table[index] + 		k%((int)(num_site_var*0.5));
        //	proc start		block row start						fast changing index
      }

      // C
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds->mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var
			+ 1 * (int)SQUARE(num_site_var/2) + k); 
        *(g.ds->mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5);
      }

      // B
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds->mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var
			+ 2 * (int)SQUARE(num_site_var/2) + k); 
        *(g.ds->mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + 0;
      }

      // D
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds->mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var
			+ 3 * (int)SQUARE(num_site_var/2) + k); 
        *(g.ds->mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5);
      }


     // computing and storing column indices Js
      if (comm_nr[dir] > 0 && op->neighbor_table[index+1+dir] >= l->num_inner_lattice_sites){

	li = boundary_table[op->neighbor_table[index + 1 + dir] % comm_nr[dir]]; //will hold the index of neighboring lattice site on neighboring domain

        for (k = 0; k < SQUARE(num_site_var/2); k ++){
//FINDME
          //A
          *(g.ds->mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) = 
			neighbors_j_start + num_site_var * li +
			 k/((int)(num_site_var*0.5));
          //C
	  *(g.ds->mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * li +
			 k/((int)(num_site_var*0.5));
          //B
          *(g.ds->mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * li +
			 k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
          //D
          *(g.ds->mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir + 1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) =
			neighbors_j_start + num_site_var * li +
			 k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
        }
	bt_index++;

        // write mu- coupling to buffer
	// also write global target and source site number to buffer
	// send both buffers
	buffer_d_pt = buffer_i_pt * num_link_var;
	//  A*
        for (k = 0; k < SQUARE(site_var/2); k ++){
          *(buff_d_send[dir] + buffer_d_pt + k) = 
				-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + k));
        }
	// -C*
        for (k = 0; k < SQUARE(site_var/2); k ++){
          *(buff_d_send[dir] + buffer_d_pt + 1 * (int)SQUARE(site_var/2) + k) = 
				 1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 1*(int)SQUARE(site_var/2) + k));
        }
	// -B*
        for (k = 0; k < SQUARE(site_var/2); k ++){
          *(buff_d_send[dir] + buffer_d_pt + 2 * (int)SQUARE(site_var/2) + k) = 
				 1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 2 * (int)SQUARE(site_var/2) + k));
        }
	// D*
        for (k = 0; k < SQUARE(site_var/2); k ++){
          *(buff_d_send[dir] + buffer_d_pt + 3 * (int)SQUARE(site_var/2) + k) = 
				-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 3 * (int)SQUARE(site_var/2) + k));
        }

	// write site index to buffer:
	*(buff_i_send[dir] + 2 * buffer_i_pt) = li;
        *(buff_i_send[dir] + 2 * buffer_i_pt + 1) = op->neighbor_table[index];
 	buffer_i_pt++;
	//printf0("node %d has neighbor in dir %d on different domain\n", node, dir);
      } else {	// neighboring lattice site is on same/my process-domain

        for (k = 0; k < SQUARE(num_site_var/2); k ++){
          //A
          *(g.ds->mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) =
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5));
          //C 
          *(g.ds->mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5));
          //B
          *(g.ds->mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
          //D
          *(g.ds->mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) =
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
        }
      }
      //printf0("hopping node: %d done!\n", node);
    }	// loop over nodes/lattice sites

  
    //printf0("hopping in dim: %d done!\n", dir);

    // sending both buffers:
    if (comm_nr[dir] > 0){

      MPI_Isend(buff_d_send[dir], comm_nr[dir] * num_link_var, MPI_COMPLEX_PRECISION,
	      l->neighbor_rank[2*dir], dir, g.comm_cart, &req);
      MPI_Isend(buff_i_send[dir], 2 * comm_nr[dir], MPI_INT, l->neighbor_rank[2*dir], dir,
	      g.comm_cart, &req);

      // HOW TO FIND NEIGHBOR?
      // int l.neighbor_rank[8] contains ranks of neighbors
      // in the order [T+ T- Z+ Z- ...]
      // no barrier here. Ordering is ensured by tag = dir
    }
    //printf0("Isends in dim: %d done!\n", dir);
  }  // loop over directions

  //printf0("hopping in mu+ done!\n");

  // mu- couplings
  for (dir = T; dir <= X; dir++){
    if (comm_nr[dir] > 0){ //there is stuff to communicate in direction dir
      MPI_Recv(buff_d_recv[dir], num_link_var * comm_nr[dir], MPI_COMPLEX_PRECISION, l->neighbor_rank[2*dir+1], dir, g.comm_cart, &s);
      MPI_Recv((buff_i_recv[dir]), 2 * comm_nr[dir], MPI_INT, l->neighbor_rank[2*dir+1], dir, g.comm_cart, &s);

      // staring column of communication-partner process
      // global rank != local rank, -> problem, since this returns global rank
      neighbors_j_start = loc_ranks[l->neighbor_rank[2*dir+1]] * l->num_inner_lattice_sites * site_var;	
      // copy buffer content to mumps_vals
      for (buffer_i_pt = 0; buffer_i_pt < comm_nr[dir]; buffer_i_pt++){
        buffer_d_pt = num_link_var * buffer_i_pt;
	// A*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(g.ds->mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + k) = 
			*(buff_d_recv[dir] + buffer_d_pt + k);
	  *(g.ds->mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5);
	  *(g.ds->mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + k) = 
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5));
        }
	// -C*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(g.ds->mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 1 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(g.ds->mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5);
	  *(g.ds->mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
        }
	// -B*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(g.ds->mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 2 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(g.ds->mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(g.ds->mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5));
        } 
	// D*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(g.ds->mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 3 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(g.ds->mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(g.ds->mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
        }
      }
    }//end if (comm_nr[dir] > 0)

    //printf0("recieved communicated couplings in dir: %d\n", dir);


    // regular mu- coupling for all nodes except communicated ones
    buffer_i_pt = 0;
    for (i = 0; i < core_end; i++){	//loop over lattice sites
      if (comm_nr[dir] > 0){
        while (i == *(buff_i_recv[dir] + 2 * buffer_i_pt + 1)){
          // skip this node because it was already communicated
          i++;
          buffer_i_pt++;
          if (i >= core_end) break;
        }
	if (i >= core_end) break;
      }
      index = 5 * i;
      // regular mu- coupling
      // A*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds->mumps_vals + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + k) = 
			-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + k));
        *(g.ds->mumps_Is   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5);
        *(g.ds->mumps_Js   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + k) = 
			(j_start + num_site_var * op->neighbor_table[index] +	 k%(int)(num_site_var*0.5));
      }
      // -C*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds->mumps_vals + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k)); 
        *(g.ds->mumps_Is   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5);
        *(g.ds->mumps_Js   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
      }
      // -B*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds->mumps_vals + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k)); 
        *(g.ds->mumps_Is   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(g.ds->mumps_Js   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5));
      }
      // D*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds->mumps_vals + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k)); 
        *(g.ds->mumps_Is   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(g.ds->mumps_Js   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
      }
    }	//loop over nodes  
    //printf0("mu- on node %d, dir %d done!\n", node, dir);
  }	//loop over directions

 
 //printf0("hopping in mu - done!\n");

#ifdef COARSE_SCALAP	//use Scalapack on coarsest level to solve.
  
//  printf0("starting coarse_scalap_setup\n"); fflush(stdout);
//  MPI_Barrier(MPI_COMM_WORLD);
  coarse_scalap_setup_PRECISION( l, threading);
//printf0("generate scalap matrix done\n");
#else
  // increase global indices by 1 to match fortran indexing.
  // spmv doesn't work then anymore
  int nnz_loc = SQUARE(site_var) * nr_nodes *9;
  for (i = 0; i < nnz_loc; i++){	//increase indices by one to match fortran indexing in MUMPS
    *(g.ds->mumps_Js + i ) = *(g.ds->mumps_Js + i ) +1;
    *(g.ds->mumps_Is + i ) = *(g.ds->mumps_Is + i ) +1;
  }
#endif


//  printf0("freeing ranks\n");
  //TODO: fix threading in here!
  //START_MASTER(threading)
  FREE( glob_ranks, int, global_comm_size);
  FREE( loc_ranks, int, global_comm_size);
//  END_MASTER(threading)
//  SYNC_CORES(threading)

//TODO: Release buffers for communication or use existing ones!

  // timing the setup
//  START_MASTER(threading)
  t1 = MPI_Wtime();
  printf0("direct coarse pre-setup time (seconds) : %f\n",t1-t0);
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  

}

void direct_solves_set_reset_PRECISION( level_struct *lx, struct Thread *threading ){

    g.coarsest_time = 0;
    g.coarsest_solve_time = 0;
    g.scalap_comm_time = 0;
    if (g.on_solve) {
      level_struct *l = lx;
      for (int i = 1; i<g.num_levels; i++) l = l->next_level;
      double t0 = 0, t1 = 0;
      if (!l->idle){ 

	  mumps_setup_PRECISION(l, threading);        //setup vals, Is, Js (and the 1D odd-even blocks for scalapack)

          START_MASTER(threading)
          t0 = MPI_Wtime();
	  END_MASTER(threading)
	  SYNC_CORES(threading)

#if defined(MUMPS_ADDS)
          printf0("starting analyze from top_level.c\n");

	  //    g.mumps_id.job = 4; //analyze and factorize
          g.mumps_id.job = 1; //analyze

          cmumps_c(&(g.mumps_id));

          printf0("analyze done, starting factorize singlethreaded from top_level.c\n");
          g.mumps_id.job = 2; //factorize
          cmumps_c(&(g.mumps_id)); //only factorize when on solve
#endif
      }
#if defined(COARSE_SCALAP)
      // collective over all processes (also idle ones): S = D_ee - D_eo D_oo^{-1} D_oe, S^{-1}
      coarse_scalap_factorize_PRECISION( l, threading );
#endif
      if (!l->idle){ 
	  t1 = MPI_Wtime();
#if defined(MUMPS_ADDS)
	  printf0("MUMPS analyze and factorize time (seconds) : %f \t from top_level.c\n",t1-t0);
#elif defined(COARSE_SCALAP)
	  printf0("Schur complement + inversion using scalapack time (seconds) : %f \t from top_level.c\n",t1-t0);
#endif
          START_MASTER(threading)
	  g.coarsest_fact_time += t1-t0;
	  END_MASTER(threading)
	  SYNC_CORES(threading)

	  //set direct solve as precond during solve phase
#if defined(MUMPS_ADDS)        
          l->p_PRECISION.preconditioner = mumps_solve_PRECISION;
	  l->p_PRECISION.eval_operator = coarse_apply_oddeven_operator_PRECISION;
	  //MUMPS uses the entire coarsest system, not only the odd-even Schur complement
	  if (g.odd_even) l->p_PRECISION.v_end *=2;
#endif
        }

#ifdef COARSE_SCALAP
      // one-time check of the direct solver (full system and Schur complement)
      coarse_scalap_selftest_PRECISION( l );
#endif
     }
}

#ifdef MUMPS_ADDS
void mumps_solve_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                            int res, level_struct *lx, struct Thread *threading )
{
    if (!lx->idle){
      START_MASTER(threading)

      g.coarsest_solve_time -= MPI_Wtime();
      //END_MASTER(threading)

      gmres_PRECISION_struct* px = &(lx->p_PRECISION);

      int i, lrank = 0;
      MPI_Comm_rank(lx->gs_PRECISION.level_comm, &lrank);


      // setting up RHS
      //TODO: if odd even is activated -> 2 * (v_end - v_start)?
      int rhs_len = lx->p_PRECISION.v_end-lx->p_PRECISION.v_start;
      for (i = 0; i < rhs_len; i++){	//set the rhs-indices to global values
	  //TODO: use correct number of processes = processes of coarsest level
        *(g.ds->mumps_irhs_loc + i) = lrank * rhs_len + i+1;		//+1 due to fortran indexing
      }

      // copying eta to local_rhs
      vector_PRECISION_copy(g.ds->mumps_rhs_loc, eta, px->v_start, px->v_end, lx );

      // centralized solution
      if (g.my_rank == 0){
	// FIXME : do some sort of casting here, to avoid warnings at compile-time
	g.mumps_id.rhs = (mumps_complex *)g.ds->mumps_SOL;
      }

      // solving
      g.mumps_id.job = 3; // solve
      cmumps_c(&(g.mumps_id));

      // distributing the solution to all processes. Must be stored in px->x
      int send_count = (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start);
      MPI_Scatter(g.ds->mumps_SOL, send_count, MPI_COMPLEX_PRECISION, phi, send_count, MPI_COMPLEX_PRECISION, 0, lx->gs_PRECISION.level_comm); // lx->gs_PRECISION.level_comm

      // counting solves and measure time not only for mumps_solve but also distributing solution to
      // processes.
      //START_MASTER(threading) 
      g.coarsest_solve_number ++;
      g.coarsest_solve_time += MPI_Wtime();
      printf0("mumps time = %f, mumps solves:  %d\n", g.coarsest_solve_time, g.coarsest_solve_number);
      END_MASTER(threading)
      SYNC_CORES(threading);
    }
}


void mumps_init_PRECISION(gmres_PRECISION_struct *p, int mumps_n, int nnz_loc, int rhs_len, level_struct *lx, Thread *threading)
{
  
    // configure MUMPS_struct
    g.mumps_id.job = JOB_INIT;
    g.mumps_id.par = 1;
    g.mumps_id.sym = 0;
    g.mumps_id.comm_fortran = (MUMPS_INT) MPI_Comm_c2f(lx->gs_PRECISION.level_comm);
    
    START_MASTER(threading)
    cmumps_c(&(g.mumps_id));
    END_MASTER(threading)
    SYNC_CORES(threading);
    
    // control parameters to define how to solve the system
    g.mumps_id.ICNTL(5) = 0;    //assembled matrix
    g.mumps_id.ICNTL(18) = 3;   //distributed local triplets for analysis and factorization
    g.mumps_id.ICNTL(20) = 10;  //distributed RHS. compare to inctl(20) = 11
    g.mumps_id.ICNTL(35) = 2;   //BLR feature is activated during factorization and solution phase
//    g.mumps_id.ICNTL(35) = 3;   //BLR feature is activablrted during factorization, not used in solve
    g.mumps_id.cntl[6] = g.mumps_drop_tol;    //dropping parameter ε    (absolute error)
    // index original/in fortran 7 but in c 6

    // linking LHS
    START_MASTER(threading)
    //printf0("setting lhs\n");
    END_MASTER(threading)
    g.mumps_id.n = mumps_n;     //needed at least on P0
    g.mumps_id.nnz_loc = nnz_loc;
    g.mumps_id.irn_loc = p->mumps_Is;
    g.mumps_id.jcn_loc = p->mumps_Js;
    g.mumps_id.a_loc = (mumps_complex *)p->mumps_vals;

    // linking RHS
    START_MASTER(threading)
    //printf0("setting rhs\n");
    END_MASTER(threading)
    g.mumps_id.nloc_rhs = rhs_len;
    g.mumps_id.rhs_loc = (mumps_complex *)p->mumps_rhs_loc;
    g.mumps_id.irhs_loc = p->mumps_irhs_loc;
    g.mumps_id.lrhs_loc = rhs_len; //leading dimension

    // solution only known to P0
    if (g.my_rank == 0){
      g.mumps_id.rhs = (mumps_complex *)p->mumps_SOL;
    }

    // control parameter for output (0 == suppressed)
    g.mumps_id.ICNTL(1) = 0;//6;        //error messages
    g.mumps_id.ICNTL(2) = 0;//1;        //diagnostic printing and statistics local to each MPI process
    g.mumps_id.ICNTL(3) = 0;//6;        //global information, collected on host (default 6)
    g.mumps_id.ICNTL(4) = 2;        //level of printing for error, warning, and diagnostic messages (default 2)

//    printf0("finished mumps_init_PRECISION\n\n\n");
    
}
#endif


#ifdef COARSE_SCALAP
/* ===========================================================================
 * Direct solver on the odd-even Schur complement of the coarsest level
 *
 *   D = [ D_ee  D_eo ]      S = D_ee - D_eo D_oo^{-1} D_oe
 *       [ D_oe  D_oo ]
 *
 * Setup (each time the coarsest operator changes):
 *   coarse_scalap_setup_PRECISION      sorts the locally assembled rows of D
 *                                      (mumps_vals/Is/Js) into D_ee, D_eo,
 *                                      D_oe and the diagonal blocks D_zz of
 *                                      the odd sites, then computes
 *                                      D_oo^{-1} D_oe locally (no communication)
 *   coarse_scalap_factorize_PRECISION  1D -> 2D block cyclic, one pgemm
 *                                      S = D_ee - D_eo (D_oo^{-1} D_oe),
 *                                      pgetrf + pgetri  ->  S^{-1} explicitly
 * Solve:
 *   coarse_scalap_solve_oddeven_PRECISION
 *       b_e <- b_e - D_eo D_oo^{-1} b_o      (sparse, existing odd-even routines)
 *       x_e  = S^{-1} b_e                    (pgemv)
 *       x_o  = D_oo^{-1} (b_o - D_oe x_e)    (sparse, existing odd-even routines)
 *
 * Every local coarsest lattice holds equally many even and odd sites, hence
 * all blocks have the global size N_h x N_h (N_h = N/2) and one descriptor per
 * layout (1D, 2D) is sufficient. In 1D every process owns the bs_h rows of its
 * own even (resp. odd) sites, in the local odd-even ordering of oe_op.
 * All memory is allocated once (direct_solver_PRECISION_setup,
 * coarse_scalap_init_PRECISION) and freed in direct_solver_PRECISION_free.
 * Only one thread per process is supported.
 * ========================================================================= */

// abort from any rank (error0 only aborts on rank 0)
static void scalap_abort_PRECISION( const char *msg, long value ) {
    printf("\x1b[31merror (rank %d): %s (%ld)\x1b[0m\n", g.my_rank, msg, value);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

// descriptor for processes which are not part of a BLACS grid
static void scalap_desc_placeholder_PRECISION( lapack_int *desc, lapack_int ictxt ) {
    desc[0] = 1;
    desc[1] = ictxt;
    for (int i = 2; i < 8; i++) desc[i] = 0;
    desc[8] = 1;
}


void coarse_scalap_init_PRECISION(level_struct *l, struct Thread *threading){

    lapack_int izero = 0, ione = 1, info = 0;
    lapack_int N_h = g.ds->N_h;      // set in direct_solver_PRECISION_setup
    lapack_int bs_h = g.ds->bs_h;

    //0. the Schur complement approach requires odd-even and equally many even
    //   and odd sites on every process of the coarsest level
    if (!g.odd_even)
	error0("the direct coarsest level solver works on the odd-even Schur complement, please enable odd-even preconditioning\n");
    if (!l->idle){
	operator_PRECISION_struct *op = &(l->oe_op_PRECISION);
	if (op->num_even_sites != op->num_odd_sites)
	    scalap_abort_PRECISION("coarse_scalap_init: number of even and odd sites differ on the coarsest level, num_even_sites =",
				   (long) op->num_even_sites);
	if ((lapack_int) op->num_even_sites * l->num_lattice_site_var != bs_h)
	    scalap_abort_PRECISION("coarse_scalap_init: inconsistent block size bs_h =", (long) bs_h);
    }

    //1. Initialize 1D - (cyclic) Blacs-Grid for reading in/out data
    lapack_int nprow = l->num_processes;
    lapack_int npcol = 1;
    lapack_int mycol;
    lapack_int ictxt1d = MPI_Comm_c2f(l->gs_PRECISION.level_comm);

    if (!l->idle){
	lapack_int *usermap = NULL;
	MALLOC(usermap, lapack_int, l->num_processes);
	for (int i = 0; i < l->num_processes; i++) usermap[i] = i;
        blacs_gridmap_( &ictxt1d, usermap, &nprow, &nprow, &npcol );
	blacs_gridinfo_( &ictxt1d, &nprow, &npcol, &g.ds->myrow, &mycol );	//set process coordinates of blacs-grid
	g.ds->blacs_ctxt1d = ictxt1d; //blacs_gridmap_ changes ictxt, therefore save new context
	FREE(usermap, lapack_int, l->num_processes);

	// every process holds the bs_h rows of its own even (resp. odd) sites
	lapack_int numr = numroc_( &N_h, &bs_h, &g.ds->myrow, &izero, &nprow );
	lapack_int lld1d = numr > 1 ? numr : 1;
	//descinit ( DESC,		 M,    N,     MB,    NB,    IRSRC, ICSRC, ICTXT,    LLD,    INFO )
	descinit_( g.ds->desc_mat1d, &N_h, &N_h,  &bs_h, &bs_h, &izero, &izero, &ictxt1d, &lld1d, &info);
	if (info != 0) error0("Error in descinit for desc_mat1d, info = %d\n", info);
	descinit_( g.ds->desc_vec1d, &N_h, &ione, &bs_h, &ione, &izero, &izero, &ictxt1d, &lld1d, &info);
	if (info != 0) error0("Error in descinit for desc_vec1d, info = %d\n", info);
    } else {
	scalap_desc_placeholder_PRECISION( g.ds->desc_mat1d, ictxt1d );
	scalap_desc_placeholder_PRECISION( g.ds->desc_vec1d, ictxt1d );
    }

    //2. Initialize 2D - cyclic Blacs-Grid for performing fast calculation
    lapack_int nprow2d = g.prow2d;
    lapack_int npcol2d = g.pcol2d;

    MPI_Comm level_comm_dup;
    MPI_Comm_dup(MPI_COMM_WORLD, &level_comm_dup);
    lapack_int ictxt2d = MPI_Comm_c2f(level_comm_dup);

    {
	lapack_int *usermap = NULL;
	MALLOC(usermap, lapack_int, g.num_processes);
	for (int i = 0; i < g.num_processes; i++) usermap[i] = i;
	blacs_gridmap_( &ictxt2d, usermap, &nprow2d, &nprow2d, &npcol2d );
	blacs_gridinfo_( &ictxt2d, &nprow2d, &npcol2d, &g.ds->myrow2d, &g.ds->mycol2d );
	FREE(usermap, lapack_int, g.num_processes);
    }
    g.ds->blacs_ctxt2d = ictxt2d; //blacs_gridmap_ changes ictxt, therefore save new context

    lapack_int bs2d = g.bs2d;
    lapack_int numr2d = 0, numc2d = 0;
    if (g.ds->myrow2d >= 0) {
	numr2d = numroc_( &N_h, &bs2d, &g.ds->myrow2d, &izero, &nprow2d );
	numc2d = numroc_( &N_h, &bs2d, &g.ds->mycol2d, &izero, &npcol2d );
	lapack_int lld2d = numr2d > 1 ? numr2d : 1;
	//descinit ( DESC,		 M,    N,     MB,    NB,    IRSRC, ICSRC, ICTXT,    LLD,    INFO )
	descinit_( g.ds->desc_mat2d, &N_h, &N_h,  &bs2d, &bs2d, &izero, &izero, &ictxt2d, &lld2d, &info);
	if (info != 0) error0("Error in descinit for desc_mat2d, info = %d\n", info);
	descinit_( g.ds->desc_vec2d, &N_h, &ione, &bs2d, &ione, &izero, &izero, &ictxt2d, &lld2d, &info);
	if (info != 0) error0("Error in descinit for desc_vec2d, info = %d\n", info);
    } else {
	scalap_desc_placeholder_PRECISION( g.ds->desc_mat2d, ictxt2d );
	scalap_desc_placeholder_PRECISION( g.ds->desc_vec2d, ictxt2d );
    }

    //3. allocate the 2D arrays once (exact local sizes from numroc)
    g.ds->len_mat2d = MAX( (long) numr2d * (long) numc2d, 1 );
    g.ds->len_vec2d = MAX( (long) numr2d, 1 );
    g.ds->len_ipiv  = (long) numr2d + bs2d;
    MALLOC( g.ds->S2d, complex_float, g.ds->len_mat2d );
    MALLOC( g.ds->Deo2d, complex_float, g.ds->len_mat2d );
    MALLOC( g.ds->DooinvDoe2d, complex_float, g.ds->len_mat2d );
    MALLOC( g.ds->rhs2d, complex_float, g.ds->len_vec2d );
    MALLOC( g.ds->sol2d, complex_float, g.ds->len_vec2d );
    MALLOC( g.ds->ipiv, lapack_int, g.ds->len_ipiv );
    memset( g.ds->S2d, 0, g.ds->len_mat2d * sizeof(complex_float) );
    memset( g.ds->Deo2d, 0, g.ds->len_mat2d * sizeof(complex_float) );
    memset( g.ds->DooinvDoe2d, 0, g.ds->len_mat2d * sizeof(complex_float) );
    memset( g.ds->rhs2d, 0, g.ds->len_vec2d * sizeof(complex_float) );
    memset( g.ds->sol2d, 0, g.ds->len_vec2d * sizeof(complex_float) );
    memset( g.ds->ipiv, 0, g.ds->len_ipiv * sizeof(lapack_int) );

    //4. workspace query for pgetri, allocate workspace once
    if (g.ds->myrow2d >= 0) {
	complex_PRECISION work_query = 0;
	lapack_int iwork_query = 0, lwork = -1, liwork = -1;
	pgetri_PRECISION( &N_h, g.ds->S2d, &ione, &ione, g.ds->desc_mat2d, g.ds->ipiv,
			  &work_query, &lwork, &iwork_query, &liwork, &info );
	if (info != 0) error0("Error during pgetri_ workspace query, info = %d\n", info);
	g.ds->lwork  = MAX( (long) creal(work_query), 1 );
	g.ds->liwork = MAX( (long) iwork_query, 1 );
	MALLOC( g.ds->work, complex_float, g.ds->lwork );
	MALLOC( g.ds->iwork, lapack_int, g.ds->liwork );
    }

    printf0("coarse_scalap_init: Schur complement of order N_h = %ld (full system N = %ld), 2D grid %ld x %ld, bs2d = %ld\n",
	    (long) N_h, 2*(long) N_h, (long) nprow2d, (long) npcol2d, (long) bs2d);
}


void coarse_scalap_setup_PRECISION(level_struct *l, struct Thread *threading){
    /* Sorts the locally assembled rows of D (mumps_vals, global row indices
     * mumps_Is, global column indices mumps_Js, local odd-even ordering of
     * oe_op: even sites first) into
     *   S1d         <- D_ee            (even rows, even columns)
     *   Deo1d       <- D_eo            (even rows, odd  columns)
     *   DooinvDoe1d <- D_oe            (odd  rows, even columns)
     *   Doo_blocks  <- D_zz            (odd  rows, odd  columns: diagonal blocks only)
     * and computes D_oo^{-1} D_oe locally. All 1D arrays are column major with
     * leading dimension bs_h, columns in the global even (resp. odd) numbering
     *   r * bs_h + (local even/odd site) * nv + k.
     */
    if (l->idle) return;

    operator_PRECISION_struct *op = &(l->oe_op_PRECISION);
    long nv = l->num_lattice_site_var, nv2 = nv*nv;
    long n = l->num_inner_lattice_sites, ne = op->num_even_sites, no = op->num_odd_sites;
    long bs = n*nv, bs_h = g.ds->bs_h, N_h = g.ds->N_h;
    long nnz_loc = 9 * nv2 * n;
    long num_misplaced = 0;
    int lrank = 0;
    MPI_Comm_rank(l->gs_PRECISION.level_comm, &lrank);

    complex_float *S1d = g.ds->S1d, *Deo1d = g.ds->Deo1d,
		  *DooinvDoe1d = g.ds->DooinvDoe1d, *Doo = g.ds->Doo_blocks;

    memset( S1d, 0, g.ds->len_mat1d * sizeof(complex_float) );
    memset( Deo1d, 0, g.ds->len_mat1d * sizeof(complex_float) );
    memset( DooinvDoe1d, 0, g.ds->len_mat1d * sizeof(complex_float) );
    memset( Doo, 0, g.ds->len_doo * sizeof(complex_float) );

    for (long e = 0; e < nnz_loc; e++){
	long r = (long) g.ds->mumps_Is[e] - lrank*bs;	// local row (same convention as before: Is = row, Js = column)
	long c = (long) g.ds->mumps_Js[e];		// global column
	if (r < 0 || r >= bs || c < 0 || c >= bs * l->num_processes) { num_misplaced++; continue; }

	long s_r = r / nv, k_r = r % nv;			// local site / component of the row
	long c_rank = c / bs, s_c = (c % bs) / nv, k_c = c % nv;	// owner, local site, component of the column
	int r_even = (s_r < ne), c_even = (s_c < ne);

	long rh = (r_even ? s_r : s_r - ne) * nv + k_r;			// local row in the even/odd half
	long ch = c_rank * bs_h + (c_even ? s_c : s_c - ne) * nv + k_c;	// global column in the even/odd numbering

	if (r_even && c_even) {			// D_ee: self coupling only
	    if (c_rank != lrank || s_c != s_r) { num_misplaced++; continue; }
	    S1d[ch*bs_h + rh] += g.ds->mumps_vals[e];
	} else if (r_even) {			// D_eo
	    Deo1d[ch*bs_h + rh] += g.ds->mumps_vals[e];
	} else if (c_even) {			// D_oe
	    DooinvDoe1d[ch*bs_h + rh] += g.ds->mumps_vals[e];
	} else {				// D_oo: self coupling only
	    if (c_rank != lrank || s_c != s_r) { num_misplaced++; continue; }
	    Doo[(s_r - ne)*nv2 + k_c*nv + k_r] += g.ds->mumps_vals[e];
	}
    }
    if (num_misplaced > 0)
	scalap_abort_PRECISION("coarse_scalap_setup: matrix entries outside the expected odd-even block structure", num_misplaced);

    // D_oo^{-1} D_oe: for every local odd site z solve D_zz X = (rows of z in D_oe), in place.
    // The rows of z form a nv x N_h submatrix of DooinvDoe1d with leading dimension bs_h.
    lapack_int ipiv_loc[nv];
    for (long z = 0; z < no; z++){
	lapack_int info = gesv_PRECISION( LAPACK_COL_MAJOR, (lapack_int) nv, (lapack_int) N_h,
					  Doo + z*nv2, (lapack_int) nv, ipiv_loc,
					  DooinvDoe1d + z*nv, (lapack_int) bs_h );
	if (info != 0) scalap_abort_PRECISION("coarse_scalap_setup: gesv for D_oo^{-1} D_oe failed, info =", (long) info);
    }
}


void coarse_scalap_factorize_PRECISION( level_struct *l, struct Thread *threading ){
    /* collective over all processes (1D and 2D grid):
     *   S2d = D_ee - D_eo (D_oo^{-1} D_oe),  S2d <- S^{-1}
     */
    lapack_int ione = 1, info = 0;
    lapack_int N_h = g.ds->N_h;

    //1. distribute the three matrices to the 2d block cyclic pattern
    scalap_1d_2d_A_PRECISION( l, threading );

    if ( g.ds->myrow2d >= 0 ){
	//2. S = D_ee - D_eo (D_oo^{-1} D_oe)
	char trans = 'N';
	complex_PRECISION alpha = -1.0, beta = 1.0;
	pgemm_PRECISION( &trans, &trans, &N_h, &N_h, &N_h, &alpha,
			 g.ds->Deo2d, &ione, &ione, g.ds->desc_mat2d,
			 g.ds->DooinvDoe2d, &ione, &ione, g.ds->desc_mat2d,
			 &beta, g.ds->S2d, &ione, &ione, g.ds->desc_mat2d );

	//3. LU decomposition of S
	pgetrf_PRECISION( &N_h, &N_h, g.ds->S2d, &ione, &ione, g.ds->desc_mat2d, g.ds->ipiv, &info );
        if (info != 0 ) error0("Error during pgetrf_(), info = %d\n", info);

	//4. explicit inverse S^{-1} (workspace allocated in coarse_scalap_init_PRECISION)
	lapack_int lwork = g.ds->lwork, liwork = g.ds->liwork;
	pgetri_PRECISION( &N_h, g.ds->S2d, &ione, &ione, g.ds->desc_mat2d, g.ds->ipiv,
			  g.ds->work, &lwork, g.ds->iwork, &liwork, &info );
        if (info != 0 ) error0("Error during pgetri_(), info = %d\n", info);
    }
}


void coarse_scalap_apply_Sinv_PRECISION( vector_PRECISION x, vector_PRECISION b,
					 level_struct *l, struct Thread *threading ){
    /* x_e = S^{-1} b_e on the even part of the DDalphaAMG vectors x and b
     * (first bs_h entries of each process). Collective over all processes.
     */
    lapack_int ione = 1;
    lapack_int N_h = g.ds->N_h;
    char trans = 'N';
    complex_PRECISION alpha = 1.0, beta = 0.0;

    //rhs2d = b_e
    scalap_1d_2d_vec_PRECISION( b, l, threading );
    //sol2d = S^{-1} rhs2d   (separate output vector, PBLAS does not allow x == y)
    if ( g.ds->myrow2d >= 0 ){
	pgemv_PRECISION( &trans, &N_h, &N_h, &alpha, g.ds->S2d, &ione, &ione, g.ds->desc_mat2d,
			 g.ds->rhs2d, &ione, &ione, g.ds->desc_vec2d, &ione,
			 &beta, g.ds->sol2d, &ione, &ione, g.ds->desc_vec2d, &ione );
    }
    //x_e = sol2d
    scalap_2d_1d_vec_PRECISION( x, l, threading );
}


void coarse_scalap_solve_oddeven_PRECISION( vector_PRECISION x, vector_PRECISION b,
					    level_struct *l, struct Thread *threading ){
    /* solves D x = b on the coarsest level via the odd-even Schur complement.
     * b is overwritten (as in coarse_solve_odd_even_PRECISION).
     * Collective over all processes, idle processes only take part in step 2.
     */
    operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

    if (!l->idle){
        START_MASTER(threading)
        g.coarsest_solve_time -= MPI_Wtime();
	END_MASTER(threading)
        SYNC_CORES(threading);

	//1. x_o = D_oo^{-1} b_o,  b_e <- b_e - D_eo x_o
	coarse_diag_oo_inv_PRECISION( x, b, op, l, threading );
	coarse_n_hopping_term_PRECISION( b, x, op, _EVEN_SITES, l, threading );
    }

    //2. x_e = S^{-1} b_e
    coarse_scalap_apply_Sinv_PRECISION( x, b, l, threading );

    if (!l->idle){
	//3. b_o <- b_o - D_oe x_e,  x_o = D_oo^{-1} b_o
	coarse_n_hopping_term_PRECISION( b, x, op, _ODD_SITES, l, threading );
	coarse_diag_oo_inv_PRECISION( x, b, op, l, threading );

        START_MASTER(threading)
        g.coarsest_solve_number ++;
        g.coarsest_solve_time += MPI_Wtime();
        printf0("scalap time  = %f, scalap solves:  %d\n", g.coarsest_solve_time, g.coarsest_solve_number);
	END_MASTER(threading)
        SYNC_CORES(threading);
    }
}


void coarse_scalap_selftest_PRECISION( level_struct *l ){
    /* executed once (collective over all processes):
     * test 1: full system      ||b - D x|| / ||b||            with x from coarse_scalap_solve_oddeven
     * test 2: Schur complement ||x_e - S^{-1} S x_e|| / ||x_e|| with S applied matrix free
     *                                                   (coarse_apply_schur_complement_PRECISION)
     */
    static int checked = 0;
    if (checked) return;
    checked = 1;

    operator_PRECISION_struct *op = &(l->oe_op_PRECISION);
    int vs = MAX(l->vector_size, 1), ivs = l->inner_vector_size;
    int nev = l->idle ? 0 : op->num_even_sites * l->num_lattice_site_var;
    vector_PRECISION b = NULL, w = NULL, x = NULL, r = NULL;
    MALLOC(b, complex_PRECISION, vs); MALLOC(w, complex_PRECISION, vs);
    MALLOC(x, complex_PRECISION, vs); MALLOC(r, complex_PRECISION, vs);
    memset(b, 0, vs*sizeof(complex_PRECISION)); memset(w, 0, vs*sizeof(complex_PRECISION));
    memset(x, 0, vs*sizeof(complex_PRECISION)); memset(r, 0, vs*sizeof(complex_PRECISION));

    // keep solve statistics untouched
    int n0 = g.coarsest_solve_number; double t0s = g.coarsest_solve_time, t0c = g.scalap_comm_time;

    //test 1: full system
    srand(1234);
    if (!l->idle) {
	vector_PRECISION_define_random(b, 0, ivs, l);
	vector_PRECISION_copy(w, b, 0, ivs, l);			// rhs is overwritten by the solve
    }
    coarse_scalap_solve_oddeven_PRECISION(x, w, l, no_threading);
    if (!l->idle) {
	coarse_diag_ee_PRECISION(r, x, op, l, no_threading);
	coarse_diag_oo_PRECISION(r, x, op, l, no_threading);
	coarse_hopping_term_PRECISION(r, x, op, _FULL_SYSTEM, l, no_threading);
	vector_PRECISION_minus(r, b, r, 0, ivs, l);		// r = b - D x
	printf0("direct solver check (full system):      ||b - D x||/||b||           = %.2e\n",
		(double)(global_norm_PRECISION(r, 0, ivs, l, no_threading) /
			 global_norm_PRECISION(b, 0, ivs, l, no_threading)));
    }

    //test 2: Schur complement, dense S^{-1} against matrix free S
    if (!l->idle) {
	vector_PRECISION_define(x, 0, 0, ivs, l);
	vector_PRECISION_define_random(x, 0, nev, l);		// random x_e, x_o = 0
	coarse_apply_schur_complement_PRECISION(w, x, op, l, no_threading);	// w_e = S x_e
    }
    coarse_scalap_apply_Sinv_PRECISION(r, w, l, no_threading);	// r_e = S^{-1} w_e
    if (!l->idle) {
	PRECISION norm_x = global_norm_PRECISION(x, 0, nev, l, no_threading);
	vector_PRECISION_minus(r, x, r, 0, nev, l);		// r_e = x_e - S^{-1} S x_e
	printf0("direct solver check (Schur complement): ||x_e - S^-1 S x_e||/||x_e|| = %.2e\n",
		(double)(global_norm_PRECISION(r, 0, nev, l, no_threading) / norm_x));
    }

    g.coarsest_solve_number = n0; g.coarsest_solve_time = t0s; g.scalap_comm_time = t0c;
    FREE(b, complex_PRECISION, vs); FREE(w, complex_PRECISION, vs);
    FREE(x, complex_PRECISION, vs); FREE(r, complex_PRECISION, vs);
}


void scalap_1d_2d_A_PRECISION(level_struct *l, struct Thread *threading){
    // copies D_ee, D_eo and D_oo^{-1} D_oe from the 1D to the 2D block cyclic pattern
    lapack_int N_h = g.ds->N_h;
    lapack_int ione = 1; //starting indices

    //pgemr2d_PRECISION( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
    if (g.ds->myrow >= 0 || g.ds->myrow2d >= 0){
        pgemr2d_PRECISION( &N_h, &N_h, g.ds->S1d, &ione, &ione, g.ds->desc_mat1d,
	    g.ds->S2d, &ione, &ione, g.ds->desc_mat2d, &(g.ds->blacs_ctxt2d) );
        pgemr2d_PRECISION( &N_h, &N_h, g.ds->Deo1d, &ione, &ione, g.ds->desc_mat1d,
	    g.ds->Deo2d, &ione, &ione, g.ds->desc_mat2d, &(g.ds->blacs_ctxt2d) );
        pgemr2d_PRECISION( &N_h, &N_h, g.ds->DooinvDoe1d, &ione, &ione, g.ds->desc_mat1d,
	    g.ds->DooinvDoe2d, &ione, &ione, g.ds->desc_mat2d, &(g.ds->blacs_ctxt2d) );
    }
}

void scalap_1d_2d_vec_PRECISION( vector_PRECISION vec, level_struct *l, struct Thread *threading){
    // copies the even part of vec (first bs_h local entries) to rhs2d

    START_MASTER(threading)
    g.scalap_comm_time -= MPI_Wtime();
    END_MASTER(threading)
    SYNC_CORES(threading);

    lapack_int N_h = g.ds->N_h;
    lapack_int ione = 1; //starting indices

    //pgemr2d_PRECISION( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
    if (g.ds->myrow >= 0 || g.ds->myrow2d >= 0){
	pgemr2d_PRECISION( &N_h, &ione, vec, &ione, &ione, g.ds->desc_vec1d,
	    g.ds->rhs2d, &ione, &ione, g.ds->desc_vec2d, &(g.ds->blacs_ctxt2d) );
    }
    START_MASTER(threading)
    g.scalap_comm_time += MPI_Wtime();
    END_MASTER(threading)
    SYNC_CORES(threading);
}

void scalap_2d_1d_vec_PRECISION( vector_PRECISION vec, level_struct *l, struct Thread *threading){
    // copies sol2d to the even part of vec (first bs_h local entries)

    START_MASTER(threading)
    g.scalap_comm_time -= MPI_Wtime();
    END_MASTER(threading)
    SYNC_CORES(threading);

    lapack_int N_h = g.ds->N_h;
    lapack_int ione = 1; //starting indices

    //pgemr2d_PRECISION( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
    if (g.ds->myrow >= 0 || g.ds->myrow2d >= 0){
	pgemr2d_PRECISION( &N_h, &ione, g.ds->sol2d, &ione, &ione, g.ds->desc_vec2d,
	    vec, &ione, &ione, g.ds->desc_vec1d, &(g.ds->blacs_ctxt2d) );
    }
    START_MASTER(threading)
    g.scalap_comm_time += MPI_Wtime();
    END_MASTER(threading)
    SYNC_CORES(threading);
}
#endif
#endif
