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

#include "mumps_PRECISION.h"

#ifdef COARSE_SCALAP
void blacs_get_(const int*, const int*, int*);
void blacs_pinfo_(int*, int*);
void blacs_gridinit_(int*, const char*, const int*, const int*);
void blacs_gridinfo_(const int*, int*, int*, int*, int*);
void descinit_(int*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, int*);
int numroc_(const int*, const int*, const int*, const int*, const int*);

void pgesv_PRECISION(const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, int*, complex_PRECISION*, const int*, const int*, const int*, int* );
void pgetrf_PRECISION(const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, int*, int* );
//		void(const int *, const int *, _Complex float *, const int *, const int *, const int *, int *, int *)
void pgetrs_PRECISION(const char*, const int*, const int*, const complex_PRECISION*, const int*, const int*, const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, int* );
//		 void(const char *, const int *, const int *, const _Complex float *, const int *, const int *, const int *, const int *, _Complex float *, const int *, const int *, const int *, int *)
void pgemv_PRECISION(char*, int*, int*, PRECISION*, PRECISION*, int*, int*, int*, PRECISION*, int*,
	int*, int*, int*, PRECISION*, PRECISION*, int*, int*, int*, int*);
#endif

void mumps_setup_PRECISION(level_struct *l, struct Thread *threading){
  /* Setting Mumps data-storage
   *
   * l will be coarsest level
   */



  SYNC_CORES(threading)
  START_MASTER(threading)
  printf0("starting setup!\n");
  
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
      *(px->mumps_Is +k) = i_start + j * site_var + (int)(i/site_var);	// col indices
      *(px->mumps_Js +k) = j_start + j * site_var + (i % site_var); 	// row indices
    }
    k += skip;
  }

  printf0("inital counting done!\n");


  // A B
  // C D    all A, B, C, D stored columnwise
  // putting values to mumps formatting for SELF COUPLING = CLOVER
  for (j = 0; j < nr_nodes; j++){
    // A
    for (k = 0, r = 0; r < num_eig_vect; r++, k++){
      for (c = 0; c < r; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] = *(clover_pt + k); //clover is
            // triangular store in DDalphaAMG
        px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] = conj_PRECISION(*(clover_pt + k));
      }
      px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] = *(clover_pt + k); // diagonal element
    }
    clover_pt += clover_step_size1; // bend pointer to next piece of memory/clover part
    // D
    for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
      for (c = num_eig_vect; c < r; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] = *(clover_pt + k);
        px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] = conj_PRECISION(*(clover_pt + k));
      }
      px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] = *(clover_pt + k);
    }
    clover_pt += clover_step_size1;
    // C
    for (r = num_eig_vect, k = 0; r < 2*num_eig_vect; r++){
      for (c = 0; c < num_eig_vect; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + (r * site_var) + c] = -1.0*(conj_PRECISION(*(clover_pt + k)));
      }
    }

    // no clover_pt correction, clover is hermitian, use same data as for C 
    // B store column-wise / transposed from former storage
    for (r = 0, k = 0; r < num_eig_vect; r++){
      for (c = 0; c < num_eig_vect; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + (c * site_var) + r + num_eig_vect] = *(clover_pt + k);
      }
    }
    clover_pt += clover_step_size2; // bend pointer to next clover for next lattice site
  }

  printf0("clover part done!\n");

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
        px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] += *(tm_block_pt + k);
        px->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] -= conj_PRECISION(*(tm_block_pt + k));
      }
      px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] += *(tm_block_pt + k);
    }
    tm_block_pt += block_step_size;

    // D
    for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
      for (c = num_eig_vect; c < r; c++, k++){
        px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] += *(tm_block_pt + k);
        px->mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] -= conj_PRECISION(*(tm_block_pt + k));
      }
      px->mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] += *(tm_block_pt + k);
    }
    tm_block_pt += block_step_size;
  }
  printf0("twisted mass part done!\n");
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

  int *glob_ranks, *loc_ranks; //will contain global ranks and corresponding local ranks

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

  printf0("allocs for hopping done!\n");

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
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir + 1)*num_link_var + k) = 
			-1.0 * *(op->D + 	num_4link_var*op->neighbor_table[index] + 	dir*num_link_var + k);
	//					find correct block row				start of mu- coupling
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) = 
			i_start +	num_site_var * op->neighbor_table[index] + 		k%((int)(num_site_var*0.5));
        //	proc start		block row start						fast changing index
      }

      // C
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k);
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5);
      }

      // B
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k);
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + 0;
      }

      // D
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k);
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5);
      }


     // computing and storing column indices Js
      if (comm_nr[dir] > 0 && op->neighbor_table[index+1+dir] >= l->num_inner_lattice_sites){

        for (k = 0; k < SQUARE(num_site_var/2); k ++){
//FINDME
          //A
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) = 
			neighbors_j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5));
          //C
	  *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5));
          //B
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
          //D
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir + 1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) =
			neighbors_j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
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
        *(buff_i_send[dir] + 2 * buffer_i_pt) = boundary_table[op->neighbor_table[index + 1 + dir] % comm_nr[dir]];
        *(buff_i_send[dir] + 2 * buffer_i_pt + 1) = op->neighbor_table[index];
 	buffer_i_pt++;
	printf0("node %d has neighbor in dir %d on different domain\n", node, dir);
      } else {	// neighboring lattice site is on same/my process-domain



        for (k = 0; k < SQUARE(num_site_var/2); k ++){
          //A
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) =
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5));
          //C 
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5));
          //B
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
          //D
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) =
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
        }
      }
      printf0("hopping node: %d done!\n", node);
    }	// loop over nodes/lattice sites

  
    printf0("hopping in dim: %d done!\n", dir);

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
    printf0("Isends in dim: %d done!\n", dir);
  }  // loop over directions






  printf0("hopping in mu+ done!\n");

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
	  *(l->p_PRECISION.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + k) = 
			*(buff_d_recv[dir] + buffer_d_pt + k);
	  *(l->p_PRECISION.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5);
	  *(l->p_PRECISION.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + k) = 
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] +
			k%(int)(num_site_var*0.5));
        }
	// -C*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_PRECISION.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 1 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(l->p_PRECISION.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5);
	  *(l->p_PRECISION.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] +
			 k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
        }
	// -B*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_PRECISION.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 2 * SQUARE((int)(num_site_var*0.5)) + k); 
	  *(l->p_PRECISION.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(l->p_PRECISION.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] +
			 k%(int)(num_site_var*0.5));
        } 
	// D*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_PRECISION.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 3 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(l->p_PRECISION.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(l->p_PRECISION.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] +
			 k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
        }
      }
    }//end if (comm_nr[dir] > 0)


    printf0("recieved communicated couplings in dir: %d\n", dir);


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
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + k) = 
			-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + k));
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5);
        *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + k) = 
			(j_start + num_site_var * op->neighbor_table[index] +
			 k%(int)(num_site_var*0.5));
      }
      // -C*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5);
        *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] +
			 k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
      }
      // -B*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] +
			 k%(int)(num_site_var*0.5));
      }
      // D*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] +
			 k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
      }
    }	//loop over nodes  
    printf0("mu- on node %d, dir %d done!\n", node, dir);
  }	//loop over directions

  
  printf0("hopping in mu - done!\n");

#ifdef COARSE_SCALAP	//use Scalapack on coarsest level to solve.
  coarse_scalap_setup_PRECISION( l, threading);
  printf0("generate scalap matrix done\n");
#else
  // increase global indices by 1 to match fortran indexing.
  // spmv doesn't work then anymore
  int nnz_loc = SQUARE(site_var) * nr_nodes *9;
  for (i = 0; i < nnz_loc; i++){	//increase indices by one to match fortran indexing in MUMPS
    *(l->p_PRECISION.mumps_Js + i ) = *(l->p_PRECISION.mumps_Js + i ) +1;
    *(l->p_PRECISION.mumps_Is + i ) = *(l->p_PRECISION.mumps_Is + i ) +1;
  }
#endif



  printf0("freeing ranks\n");
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
        *(px->mumps_irhs_loc + i) = lrank * rhs_len + i+1;		//+1 due to fortran indexing
      }

      // copying eta to local_rhs
      vector_PRECISION_copy(px->mumps_rhs_loc, eta, px->v_start, px->v_end, lx );

      // centralized solution
      if (g.my_rank == 0){
	// FIXME : do some sort of casting here, to avoid warnings at compile-time
	g.mumps_id.rhs = px->mumps_SOL;
      }

      // solving
      g.mumps_id.job = 3; // solve
      cmumps_c(&(g.mumps_id));

      // distributing the solution to all processes. Must be stored in px->x
      int send_count = (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start);
      MPI_Scatter(px->mumps_SOL, send_count, MPI_COMPLEX_PRECISION, phi, send_count, MPI_COMPLEX_PRECISION, 0, lx->gs_PRECISION.level_comm); // lx->gs_PRECISION.level_comm

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
    printf0("setting lhs\n");
    END_MASTER(threading)
    g.mumps_id.n = mumps_n;     //needed at least on P0
    g.mumps_id.nnz_loc = nnz_loc;
    g.mumps_id.irn_loc = p->mumps_Is;
    g.mumps_id.jcn_loc = p->mumps_Js;
    g.mumps_id.a_loc = p->mumps_vals;

    // linking RHS
    START_MASTER(threading)
    printf0("setting rhs\n");
    END_MASTER(threading)
    g.mumps_id.nloc_rhs = rhs_len;
    g.mumps_id.rhs_loc = p->mumps_rhs_loc;
    g.mumps_id.irhs_loc = p->mumps_irhs_loc;
    g.mumps_id.lrhs_loc = rhs_len; //leading dimension

    // solution only known to P0
    if (g.my_rank == 0){
      g.mumps_id.rhs = p->mumps_SOL;
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
void coarse_scalap_solve_PRECISION(vector_PRECISION phi, vector_PRECISION Dphi,
                           vector_PRECISION eta, int res, level_struct *l,
                           struct Thread *threading){
    if (!l->idle){



	printf0("T %d, Z  %d, Y %d, X%d\n", T, Z, Y, X);
	for (int i = 0; i< g.num_processes; i++){
	    if(g.my_rank == i){
		printf("my rank: %d, my coords: 0:%d 1:%d 2:%d 3:%d \t X:%d, Y:%d, Z:%d, T:%d\n", g.my_rank,
		    g.my_coords[0], g.my_coords[1], g.my_coords[2], g.my_coords[3], 
		    g.my_coords[X], g.my_coords[Y], g.my_coords[Z], g.my_coords[T] ); 
		fflush(stdout);
	    }	    
	    MPI_Barrier(MPI_COMM_WORLD);
	}

	printf0("\n\n\n");
	printf0("local_lattice: 0 %d, 1 %d, 2 %d, 3 %d, \t X %d, Y %d, Z %d, T %d\n",
		l->local_lattice[0], l->local_lattice[1], l->local_lattice[2], l->local_lattice[3],
		l->local_lattice[X], l->local_lattice[Y], l->local_lattice[Z], l->local_lattice[T]);

	MPI_Barrier(MPI_COMM_WORLD);




	vector_PRECISION test2 = NULL; 
	MALLOC(test2, complex_PRECISION, l->inner_vector_size);
	memset(test2, 0, l->inner_vector_size * sizeof(complex_PRECISION));

	int rax = 2, ps = 5;
	if (g.my_rank != rax) memset(eta, 0, l->inner_vector_size * sizeof(complex_PRECISION));

	vector_PRECISION_copy(test2, eta, 0, l->inner_vector_size, l);
	
	
	
	for (int ra = 0; ra < ps; ra++){
	    if (g.my_rank == ra) for (int i = 0; i < l->num_inner_lattice_sites * l->num_lattice_site_var; i++){
		if( creal(eta[i]) == 0) printf("P%d, i: %5d, site: %5d\t:\n", g.my_rank, i, i/l->num_lattice_site_var);		
		else printf("P%d, i: %5d, site: %5d\t: %+f %+fi\n", g.my_rank, i, i/l->num_lattice_site_var, CSPLIT(eta[i]));		
		i = i + l->num_lattice_site_var - 1;
	    }
	    fflush(stdout);
	    MPI_Barrier(MPI_COMM_WORLD);
	}

	
        translate2scalap_vectors_PRECISION( l, eta);
	printf0("hinweg abgeschlossen!\n");

	for (int ra = 0; ra < ps; ra++){
	    if (g.my_rank == ra) for (int i = 0; i < l->num_inner_lattice_sites * l->num_lattice_site_var; i++){
		if ( creal(eta[i]) == 0) printf("P%d, i: %5d, site: %5d\t:\n", g.my_rank, i, i/l->num_lattice_site_var);		
		else printf("P%d, i: %5d, site: %5d\t: %+f %+fi\n", g.my_rank, i, i/l->num_lattice_site_var, CSPLIT(eta[i]));		
		i = i + l->num_lattice_site_var - 1;
	    }
	    fflush(stdout);
	    MPI_Barrier(MPI_COMM_WORLD);
	}

	
	translate2original_vectors_PRECISION( l, eta);

	printf0("rückweg abgeschlossen!\n");

	for (int ra = 0; ra < ps; ra++){
	    if (g.my_rank == ra) for (int i = 0; i < l->num_inner_lattice_sites * l->num_lattice_site_var; i++){
		if ( creal(eta[i]) == 0) printf("P%d, i: %5d, site: %5d\t:\n", g.my_rank, i, i/l->num_lattice_site_var);		
		else printf("P%d, i: %5d, site: %5d\t: %+f %+fi\n", g.my_rank, i, i/l->num_lattice_site_var, CSPLIT(eta[i]));		
		i = i + l->num_lattice_site_var - 1;
	    }
	    fflush(stdout);
	    MPI_Barrier(MPI_COMM_WORLD);
	}


	vector_PRECISION_minus( eta, eta, test2, 0, l->inner_vector_size, l );

	printf0("DIFF:\n");

	for (int ra = 0; ra < ps; ra++){
	    if (g.my_rank == ra) for (int i = 0; i < l->num_inner_lattice_sites * l->num_lattice_site_var; i++){
		if ( creal(eta[i]) == 0) printf("P%d, i: %5d, site: %5d\t:\n", g.my_rank, i, i/l->num_lattice_site_var);		
		else printf("P%d, i: %5d, site: %5d\t: %+f %+fi\n", g.my_rank, i, i/l->num_lattice_site_var, CSPLIT(eta[i]));		
		i = i + l->num_lattice_site_var - 1;
	    }
	    fflush(stdout);
	    MPI_Barrier(MPI_COMM_WORLD);
	}



	PRECISION r2 = global_norm_PRECISION( test2, 0, l->inner_vector_size, l, threading );
	PRECISION r = global_norm_PRECISION( eta, 0, l->inner_vector_size, l, threading);
	printf0("global norm: %e\n", r/r2);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
        exit(0);

   
//	printf0("solving scalap...\n");
        START_MASTER(threading)
        g.coarsest_solve_time -= MPI_Wtime();
	END_MASTER(threading)
        SYNC_CORES(threading);
	
	int ione = 1, info = 0;
	char trans = 'N';

	int N = l->num_processes * l->num_inner_lattice_sites * l->num_lattice_site_var;
	vector_PRECISION test; 
	MALLOC(test, complex_PRECISION, l->inner_vector_size);
	memset(test, 0, l->inner_vector_size * sizeof(complex_PRECISION));
	vector_PRECISION_copy(test, eta, 0, l->inner_vector_size, l);

    //void pgetrs_PRECISION( TRANS, N, NRHS, A,	     IA, JA,	 DESCA, IPIV, B, IB, JB, DESCB, INFO );
	pgetrs_PRECISION( &trans, &N, &ione, l->p_PRECISION.dense_vals, &ione, &ione,
		l->p_PRECISION.desc_dense_vals, l->p_PRECISION.ipiv,
		test, &ione, &ione, l->p_PRECISION.desc_rhs, &info );

	// test <- LU * test

        START_MASTER(threading)
	if (info != 0 ) error0("Error during pgetrs_(), info = %d\n", info);

	//vector_copy eta -> phi 
//uncomment this line to use solution        vector_PRECISION_copy(phi, test, l->p_PRECISION.v_start, l->p_PRECISION.v_end, l );



	apply_coarse_operator_PRECISION(phi, test, l->p_PRECISION.op, l, threading);
	vector_PRECISION_minus( phi, phi, eta, 0, l->inner_vector_size, l );
	PRECISION r2 = global_norm_PRECISION( eta, 0, l->inner_vector_size, l, threading );
	PRECISION r = global_norm_PRECISION(phi, 0, l->inner_vector_size, l, threading);
	printf0("global norm: %e\n", r/r2);
	MPI_Barrier(MPI_COMM_WORLD);
	exit(0);


        g.coarsest_solve_number ++;
        g.coarsest_solve_time += MPI_Wtime();
        printf0("scalap time  = %f, scalap solves:  %d\n", g.coarsest_solve_time, g.coarsest_solve_number);
	FREE(test, complex_PRECISION, l->inner_vector_size);
	END_MASTER(threading)
        SYNC_CORES(threading);
    }
}


void coarse_scalap_factorize_PRECISION( level_struct *l, vector_PRECISION A, int* descA, int* ipiv, struct Thread *threading){
    printf0("starting scalap routine\n");

//    int izero = 0, ione = 1;
    int info = 0;


    int ia = 1, ja = 1; //starting indices (global)
    int ib = 1, jb = 1;
    int N = l->num_inner_lattice_sites * l->num_lattice_site_var * l->num_processes; 
    //		    ( M,    N,	    A,		IA, JA, DESCA, IPIV, INFO )
//TODO: Re-enable this line!    pgetrf_PRECISION( &N, &N, A, &ia, &ja, descA, ipiv, &info );    
    if (info != 0 ) error0("Error during pgetrf_(), info = %d\n", info);
    
    /*
    char trans = 'N';
    PRECISION alpha = 1.0, beta = 0.0;
    vector_PRECISION outvector;
    MALLOC( outvector, complex_PRECISION, l->num_inner_lattice_sites * l->num_lattice_site_var);
    memset( outvector, 0, l->num_inner_lattice_sites * l->num_lattice_site_var * sizeof(complex_PRECISION));


//    pgemv_PRECISION(char*, int*, int*, PRECISION*, PRECISION*, int*, int*, int*, PRECISION*, int*, int*, int*, int*, PRECISION*, PRECISION*, int*, int*, int*, int*);
    pgemv_PRECISION( &trans, &N, &N, &alpha, A, &ia, &ja, descA, B, &ib, &jb, descB, &ione, &beta,
	    outvector, &ib, &jb, descB, &ione);
    vector_PRECISION_copy(B, outvector, 0, l->num_inner_lattice_sites * l->num_lattice_site_var, l);
*/
    /*
    printf0("calling pdgesv_() .... \n");
    pgesv_PRECISION( &N, &ione, A, &ia, &ja, descA, ipiv, B, &ib, &jb, descB, &info);
    printf0("pdgesv_() done \n");
    if (info != 0 ) error0("Error during pdgesv_(), info = %d\n", info);
    */
    printf0("scalap routine done!\n");
}


void coarse_scalap_setup_PRECISION(level_struct *l, struct Thread *threading){

    //TODO: the following implementation works only for one process. For multiprocessing check r and
    //c being global and use a smart modulo operation

    int r, c; //for row and column index for a given matrix element
    memset(l->p_PRECISION.dense_vals, 0, l->num_inner_lattice_sites * l->num_lattice_site_var *
	    l->num_inner_lattice_sites * l->num_lattice_site_var * l->num_processes *
	    sizeof(complex_PRECISION));

    for (int llsite = 0; llsite < l->num_inner_lattice_sites; llsite++ ){ //loop over lattice sites
    //each lattice site contains 9 blocks of each SQUARE(site_var) elements, 1 for self coupling, +
    //2x4 hopping terms (2 in each dimension)
	for (int d = 0; d < 9; d++){ //loop over each block as mentioned above
	    for (int i = 0; i < SQUARE(l->num_lattice_site_var); i++){ //local index to copy elementwise within a block
	        r = l->p_PRECISION.mumps_Is[llsite*9*SQUARE(l->num_lattice_site_var) + d * SQUARE(l->num_lattice_site_var) + i] -
			g.my_rank * l->num_inner_lattice_sites * l->num_lattice_site_var;
	        c = l->p_PRECISION.mumps_Js[llsite*9*SQUARE(l->num_lattice_site_var) + d * SQUARE(l->num_lattice_site_var) + i];
		l->p_PRECISION.dense_vals[c * l->num_inner_lattice_sites * l->num_lattice_site_var + r] += l->p_PRECISION.mumps_vals[llsite*9*SQUARE(l->num_lattice_site_var) + d * SQUARE(l->num_lattice_site_var) + i];
	    }
	    r = l->p_PRECISION.mumps_Is[llsite*9*SQUARE(l->num_lattice_site_var) + d * SQUARE(l->num_lattice_site_var)] -
			g.my_rank * l->num_inner_lattice_sites * l->num_lattice_site_var;
	    c = l->p_PRECISION.mumps_Js[llsite*9*SQUARE(l->num_lattice_site_var) + d * SQUARE(l->num_lattice_site_var)];
	}
    }
    //REMEMBER: Column major storage for Scalapack due to fortran calls!
}


void coarse_scalap_init_PRECISION(level_struct *l, struct Thread *threading){

    int izero = 0;
    int ione = 1;
    int nprow = l->num_processes;
    int npcol = 1;
    char layout = 'R';
    int info = 0;
    int iam = 0, nprocs = 0;
    int ictxt = l->p_PRECISION.blacs_ctxt, myrow, mycol;
    int N = l->num_processes * l->num_inner_lattice_sites * l->num_lattice_site_var;

    blacs_pinfo_( &iam, &nprocs); //setting rank and number of blacs-processes
    blacs_get_(&izero, &izero, &ictxt);			//create context
    blacs_gridinit_(&ictxt, &layout, &nprow, &npcol );	//create blacs grid
    blacs_gridinfo_(&ictxt, &nprow, &npcol, &myrow, &mycol);	//set process coordinates of blacs-grid

//setting the descriptors:
	//TODO may choose a different blocksize, but in beginning start with bs = 2 * num. testvecs
	//* num inner lattice sites,
	//probably a good choice

    int bs = l->num_lattice_site_var * l->num_inner_lattice_sites;
    int nrhs = 1;
	   
    int numr = numroc_( &N, &bs, &iam, &izero, &nprow ); // number of rows stored in each process
    int lddA = numr > 1? numr : 1;	//leading dimension in A (remember, matrix elements are	stored in a column major order)

    printf0("setting the descriptor of coarsest matrix for scalapack.\n");
    //descinit ( DESC,			    M, N,    MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO )
    descinit_( l->p_PRECISION.desc_dense_vals, &N, &N, &bs, &bs, &izero, &izero, &ictxt, &lddA, &info);
    printf0("matrix descriptor done.\n");
    if (info != 0) error0("Error in descinit for DescA, info = %d\n", info);

    printf0("setting the descriptor of RHS for scalapack.\n");
    //descinit ( DESC,		    M, N,	MB, NB,	    IRSRC, ICSRC, ICTXT, LLD, INFO )
    descinit_( l->p_PRECISION.desc_rhs, &N, &nrhs, &bs, &ione, &izero, &izero, &ictxt, &lddA, &info);
    printf0("RHS descriptor done.\n");
    if (info != 0) error0("Error in descinit for DescB, info = %d\n", info);
}

void translate2scalap_vectors_PRECISION( level_struct *l, vector_PRECISION phi){
/**
  * Each process checks its own local lattice domain. 
  * For each lattice site, we find the global index
    -> j coords: process-coord in process grid * size of local domain in dir j + j coord of site
  * with the global index of the lattice site, we see the position of this "chunk" of the vector in
    -> int lex_index( int t, int z, int y, int x, int N[3] ) gives the lexicographic index
    -> check for each site whether this is in "my part" of the global vector, i.e. compare with
    something proportional to g.my_rank
  * the global "block row format" storage scheme. 
    -> if not in "my part" send it to the process whose part it is.
    -> if yes, put it in correct position, by finding the correct position i.e. the local
    lexicographic index. (use coords from outer for loop)
  * check vector, compare from where to take data, if not own domain, receive data from corresponding process (not necessary neighbour process!).
    
  */


    /* mu = [X, Y, Z, T];
    size of global lattice on this level: l->global_lattice[mu];
    size of local domain: l->local_lattice[mu];
    quota global/local: l->splitting[mu];
    number of local lattice sites: l->num_inner_lattice_sites
    global rank: g.my_rank
    local process coords: g.my_coords[mu]

       */

    
    
    vector_PRECISION phi_out = NULL; 
    MALLOC(phi_out, complex_PRECISION, l->inner_vector_size);
    memset(phi_out, 0, l->inner_vector_size * sizeof(complex_PRECISION));

    int* visited_sites = NULL;
    MALLOC(visited_sites, int, l->inner_vector_size);
    memset(visited_sites, 0, l->inner_vector_size * sizeof(int));
    complex_PRECISION* buf = NULL;
    MALLOC(buf, complex_PRECISION, l->num_lattice_site_var);
    memset(buf, 0, l->num_lattice_site_var * sizeof(complex_PRECISION));


    int ra = 0;
    
    
    MPI_Request req;
    MPI_Status s;
    int source, target; //process ranks for communication (message will be sent from source to target)
    int gx, gy, gz, gt; //site coordinates in a global scheme
    int glex, llex; //lexicographic index for site with coords(gx, gy, gz, gt)


    //process coords are in order: Z Y X T
    //local lattice is in order: X Y Z T

    
    //TODO: remove this line:
    if (g.my_rank == 0) 


    
    for (int lt = 0; lt < l->local_lattice[T]; lt++){
	gt = g.my_coords[T] * l->local_lattice[T] + lt;
	for (int lz = 0; lz < l->local_lattice[Z]; lz++){
	    gz = g.my_coords[Z] * l->local_lattice[Z] + lz;
	    for (int ly = 0; ly < l->local_lattice[Y]; ly++){
		gy = g.my_coords[Y] * l->local_lattice[Y] + ly;
		for (int lx = 0; lx < l->local_lattice[X]; lx++){
		    gx = g.my_coords[X] * l->local_lattice[X] + lx;
		   
		    glex = lex_mod_index( gt, gz, gy, gx, l->global_lattice );
	     
		    llex = lex_mod_index( lt, lz, ly, lx, l->local_lattice );

		    if (glex/l->num_inner_lattice_sites == g.my_rank ){ //my chunk of RHS
			// put in pos: glex%l->nummer_inner_lattice_sites
			for (int j = 0; j < l->num_lattice_site_var; j++){
			    phi_out[ glex%l->num_inner_lattice_sites * l->num_lattice_site_var + j]
				= phi[ llex * l->num_lattice_site_var + j];
			}
			
			visited_sites[glex%l->num_inner_lattice_sites] = 1;
			if (g.my_rank == ra) {
			    printf("   ");
			}
		    } else {
			if (g.my_rank == ra) {
			    printf("com");
			}	
			//communicate to process with rank
			target = glex/l->num_inner_lattice_sites;
			//copy corresponding data	

			//TODO: keep this line?
			memset(buf, 0, l->num_lattice_site_var * sizeof(complex_PRECISION));


			vector_PRECISION_copy( buf, phi + llex * l->num_lattice_site_var, 0, l->num_lattice_site_var, l);
			//send data (non blocking)
			MPI_Isend(buf, l->num_lattice_site_var, MPI_COMPLEX_PRECISION, target, glex, g.comm_cart, &req);
			
		    }
		    if (g.my_rank == ra) {
			printf("rank: %d, gx: %d, gy: %d, gz: %d, gt: %d, glex/gsite %d, llex %d, l/ns: %d, lx: %d, ly: %d, lz: %d, lt: %d\n",
					g.my_rank, gx, gy, gz, gt, glex, llex, glex/l->num_inner_lattice_sites,
					lx, ly, lz, lt );
		    }

		}
	    }
	}
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf0("sending done!\n");

    //check my chunk of RHS in search for open position, which are not yet
    //set/communicated. 
    int gsite; //global site index
    int coords[4];  //coordinates of lattice site
    int pcoords[4]; //coordinates of processor
   


    ra = 2;

    coords[0] = 0; coords[1] = 0; coords[2] = 0; coords[3] = 0;

//    printf0("num sites: %d\nX: %d, Y: %d, Z: %d, T:%d\n", l->num_inner_lattice_sites, X, Y, Z, T);
//    printf0("local_lattice: %d, %d, %d, %d\n", l->local_lattice[X], l->local_lattice[Y], l->local_lattice[Z], l->local_lattice[T]);
    for (int lsite = 0; lsite < l->num_inner_lattice_sites; lsite++){
	if (visited_sites[lsite] == 0){	//not yet set, receive this from other process
	    //find global site index, get coordinates from there
	    // -> get process coordinates from site coordinates
	    // -> get process rank from process coordinates
	    
	    //global index for lattice site
	    gsite = g.my_rank * l->num_inner_lattice_sites + lsite;
	    //global coordinates of lattice site
	    
	    
	    coords4d(coords, gsite, l->global_lattice);
//	    llex = lex_index( coords[T], coords[Z], coords[Y], coords[X], l->local_lattice );
//	    llex = coords[X] + l->local_lattice[X] * (coords[Y] + l->local_lattice[Y] * (coords[Z] + l->local_lattice[Z] * coords[T]));
	    llex = lsite; 
	    if (g.my_rank == ra) {
		printf("site coords: %d %d %d %d, gsite %3d, llex %3d,\t ll %d %d %d %d\n", coords[X], coords[Y], coords[Z],
			coords[T], gsite, llex, l->local_lattice[X], l->local_lattice[Y], l->local_lattice[Z], l->local_lattice[T]);
	    }
	    //lexicographic index will be used as Tag during communication
//	    glex = lex_index( coords[T], coords[Z], coords[Y], coords[X], l->global_lattice );

	    //get from site coords the process coords:
	    pcoords[T] = coords[T] / l->local_lattice[T];
	    pcoords[Z] = coords[Z] / l->local_lattice[Z];
	    pcoords[Y] = coords[Y] / l->local_lattice[Y];
	    pcoords[X] = coords[X] / l->local_lattice[X];

	    //what is the rank of process with pcoords? -> store in source
	    MPI_Cart_rank(g.comm_cart, pcoords, &source);
	    //must be different from own rank!


	    //tag = gsite (to ensure, the correct lattice site is received and there is no overtaking by other messages
	    if (g.my_rank == ra) {
	/*	printf("ll: %d, %d, %d, %d, coords: %d, %d, %d, %d, pcoords: %d, %d, %d, %d,\n",
			l->local_lattice[X], l->local_lattice[Y], l->local_lattice[Z], l->local_lattice[T],
			coords[X], coords[Y], coords[Z], coords[T], 
			pcoords[X], pcoords[Y], pcoords[Z], pcoords[T]); */
		printf("receiv. glex/gsite %d, from process %d at %d,  ... lsite %d,", gsite, source,
			g.my_rank, lsite);fflush(stdout);
		    }

	
	    
	    //TODO: keep this line?
	    memset(buf, 0, l->num_lattice_site_var * sizeof(complex_PRECISION));


    
	    //TODO: remove this line:
	    if (source == 0) 

	    MPI_Recv( buf, l->num_lattice_site_var, MPI_COMPLEX_PRECISION, source, gsite, g.comm_cart, &s);
	    //copy received data to output
	    vector_PRECISION_copy( phi_out + llex*l->num_lattice_site_var, buf, 0, l->num_lattice_site_var, l);
	
	    if (g.my_rank == ra) printf("received!\n");
	    
 

	    visited_sites[lsite] = 1;
	} //else: already set since this was in own domain.
    }


    MPI_Barrier(MPI_COMM_WORLD);
    //copy phi_out to phi
    vector_PRECISION_copy( phi, phi_out, 0, l->num_lattice_site_var * l->num_inner_lattice_sites, l);

    FREE(phi_out, complex_PRECISION, l->inner_vector_size);
    FREE(visited_sites, int, l->inner_vector_size);
    FREE(buf, complex_PRECISION, l->num_lattice_site_var);

}

void translate2original_vectors_PRECISION( level_struct *l, vector_PRECISION phi){
/**
  * This function returns from Block-Row-Storage (used in SCALAPACK) back to DDalphaAMG storage
  * scheme.

  * Each process goes over local piece of vector, check, if site belongs to own domain
  * For each lattice site, we find the global index
    -> if not in "my part" send it to the process whose part it is.
    -> if yes, put it in correct position, by finding the correct position i.e. the local
    lexicographic index. (use coords from outer for loop)
  * check own domain for not yet set lattice sites, check coordinates to know where to take data
  * from, if not own domain, receive data from corresponding process (not necessary neighbour process!).

  */


    /* mu = [X, Y, Z, T];
    size of global lattice on this level: l->global_lattice[mu];
    size of local domain: l->local_lattice[mu];
    quota global/local: l->splitting[mu];
    number of local lattice sites: l->num_inner_lattice_sites
    global rank: g.my_rank
    local process coords: g.my_coords[mu]

       */
    vector_PRECISION phi_out = NULL; 
    MALLOC(phi_out, complex_PRECISION, l->inner_vector_size);
    memset(phi_out, 0, l->inner_vector_size * sizeof(complex_PRECISION));

    int* visited_sites = NULL;
    MALLOC(visited_sites, int, l->inner_vector_size);
    memset(visited_sites, 0, l->inner_vector_size * sizeof(int));
    complex_PRECISION* buf = NULL;
    MALLOC(buf, complex_PRECISION, l->num_lattice_site_var);
    memset(buf, 0, l->num_lattice_site_var * sizeof(complex_PRECISION));

    
    MPI_Request req;
    MPI_Status s;
    int source, target; //process ranks for communication (message will be sent from source to target)
    int gx, gy, gz, gt; //site coordinates in a global scheme
    int glex; //global lexicographic index for site with coords(gx, gy, gz, gt)
    int llex; //local lexicographic index for site with coords(lx, ly, lz, lt)
    int coords[4];  //global coordinates of lattice site
    int lcoords[4];  //local coordinates of lattice site
    int pcoords[4]; //coordinates of processor
    

    coords[0] = 0; coords[1] = 0; coords[2] = 0; coords[3] = 0;
    lcoords[0] = 0; lcoords[1] = 0; lcoords[2] = 0; lcoords[3] = 0;




    int ra = 0; 



    for (int lsite = 0; lsite < l->num_inner_lattice_sites; lsite++){
	//find global site index, get coordinates from there
	// -> get process coordinates from site coordinates
	// -> get process rank from process coordinates
	
	//global index for lattice site
	glex = g.my_rank * l->num_inner_lattice_sites + lsite;
	//global coordinates of lattice site
	coords4d(coords, glex, l->global_lattice);

	//get from site coords the target process coords:
	pcoords[T] = coords[T] / l->local_lattice[T];
	pcoords[Z] = coords[Z] / l->local_lattice[Z];
	pcoords[Y] = coords[Y] / l->local_lattice[Y];
	pcoords[X] = coords[X] / l->local_lattice[X];

	//what is the rank of process with pcoords? -> store in target
	MPI_Cart_rank(g.comm_cart, pcoords, &target);

	if (g.my_rank == ra) {
	    printf("r %d, ", g.my_rank);
	}
    
	
	lcoords[T] = coords[T]%l->local_lattice[T];
	lcoords[Z] = coords[Z]%l->local_lattice[Z];
	lcoords[Y] = coords[Y]%l->local_lattice[Y];
	lcoords[X] = coords[X]%l->local_lattice[X];


	//llex = lex_mod_index(coords[T], coords[Z], coords[Y], coords[X], l->local_lattice);
	llex = lcoords[X] + l->local_lattice[X]*(lcoords[Y] + l->local_lattice[Y]*(lcoords[Z] + l->local_lattice[Z]*lcoords[T]));


	if (target == g.my_rank){	//site belongs to "my domain"
		//local coordinates 
	    lcoords[T] = coords[T]%l->local_lattice[T];
	    lcoords[Z] = coords[Z]%l->local_lattice[Z];
	    lcoords[Y] = coords[Y]%l->local_lattice[Y];
	    lcoords[X] = coords[X]%l->local_lattice[X];
	    //local lex index:
	    llex = lcoords[X] + l->local_lattice[X]*(lcoords[Y] + l->local_lattice[Y]*(lcoords[Z] + l->local_lattice[Z]*lcoords[T]));
    //this must be the same as lsite  = lex_mod_index( coords[3], coords[2], coords[1], coords[0], l->local_lattice );
	    for (int j = 0; j < l->num_lattice_site_var; j++){
		phi_out[ llex * l->num_lattice_site_var + j] = phi[ lsite * l->num_lattice_site_var + j];
	    }
	    visited_sites[ llex ] = 1;

	    if (g.my_rank == ra){
		printf("___ glex %3d, from %d, to  , lsite: %d, vs: %d, llex %d, gcoords %d %d %d %d, \t lcoords %d %d %d %d, \tpcoords %d %d %d %d\n", 
			glex, g.my_rank, lsite, visited_sites[lsite], llex, 
			coords[X], coords[Y], coords[Z], coords[T], 
			lcoords[X], lcoords[Y], lcoords[Z], lcoords[T], 
			pcoords[X], pcoords[Y], pcoords[Z], pcoords[T]);
	    }
	} else { //site belongs to domain of process with rank _target_
	    //TODO: keep this line?
	    memset(buf, 0, l->num_lattice_site_var * sizeof(complex_PRECISION));


	    //copy corresponding data	
	    if (g.my_rank == ra && target == 0 && creal(phi[lsite*l->num_lattice_site_var]) == 0) {printf("copying zero!\n"); exit(0);}
	    vector_PRECISION_copy( buf, phi + lsite * l->num_lattice_site_var, 0, l->num_lattice_site_var, l);
	    //send data (non blocking)
	    MPI_Isend(buf, l->num_lattice_site_var, MPI_COMPLEX_PRECISION, target, glex, g.comm_cart, &req);
	    if (g.my_rank == ra){
		printf("com glex %3d, from %d, to %d, lsite: %d, vs: %d, llex %d, gcoords %d %d %d %d, \t lcoords %d %d %d %d, \tpcoords %d %d %d %d\n", 
			glex, g.my_rank, target, lsite, visited_sites[lsite], llex, 
			coords[X], coords[Y], coords[Z], coords[T], 
			lcoords[X], lcoords[Y], lcoords[Z], lcoords[T], 
			pcoords[X], pcoords[Y], pcoords[Z], pcoords[T]);
	    }

	}
    }
   
    ra = 0;

    MPI_Barrier(MPI_COMM_WORLD);
    printf0("sending done!\n");

    //process coords are in order: Z Y X T
    //local lattice is in order: X Y Z T
    int lsite = 0;
    for (int lt = 0; lt < l->local_lattice[T]; lt++){
	gt = g.my_coords[T] * l->local_lattice[T] + lt;
	for (int lz = 0; lz < l->local_lattice[Z]; lz++){
	    gz = g.my_coords[Z] * l->local_lattice[Z] + lz;
	    for (int ly = 0; ly < l->local_lattice[Y]; ly++){
		gy = g.my_coords[Y] * l->local_lattice[Y] + ly;
		for (int lx = 0; lx < l->local_lattice[X]; lx++){
		    gx = g.my_coords[X] * l->local_lattice[X] + lx;

		    lsite = lex_mod_index(lt, lz, ly, lx, l->local_lattice);
//		    lsite = lx + l->local_lattice[X]*(ly + l->local_lattice[Y]*(lz + l->local_lattice[Z]*lt));
		    
		    glex = lex_mod_index(gt, gz, gy, gx, l->global_lattice);
//		    glex = gx + l->global_lattice[X]*(gy + l->global_lattice[Y]*(gz + l->global_lattice[Z] * gt)); //lex_mod_index( gt, gz, gy, gx, l->global_lattice );

		    
		    if (g.my_rank == ra) {
			    printf("rank %d, lsite %3d, glex %3d, glex_mod %3d, gt %d, gz %d, gy %d, gx %d, \t\t px %d, py %d, pz %d, pt %d\n", g.my_rank, lsite, glex, glex%l->num_inner_lattice_sites, gt, gz, gy, gx, g.my_coords[X], g.my_coords[Y], g.my_coords[Z], g.my_coords[T]);
//			    printf(" product: %d, lz: %d\n", g.my_coords[0] * l->local_lattice[2], lz);
		    }


		    if (visited_sites[ lsite ] == 0){ //not yet set, data comes from different process
/*
			// get process coords from site coords
			//process coords are in order: Z Y X T
			//local lattice is in order: X Y Z T
			pcoords[2] = gx / l->local_lattice[0];
			pcoords[1] = gy / l->local_lattice[1];
			pcoords[0] = gz / l->local_lattice[2];
			pcoords[3] = gt / l->local_lattice[3];

			
			//what is the rank of process with pcoords? -> store in source
			MPI_Cart_rank(g.comm_cart, pcoords, &source);  //must be different from my_rank
*/

			source = glex/l->num_inner_lattice_sites;

			if (g.my_rank == ra) {
//				printf("pX %d, pY %d, pZ %d, pT %d, gx %d, llx %d, \tgy %d lly %d, \tgz %d, llz %d, \tgt %d, llt %d \t", pcoords[2], pcoords[1], pcoords[0], pcoords[3], gx, l->local_lattice[0], gy, l->local_lattice[1], gz, l->local_lattice[2], gt, l->local_lattice[3]);
				printf("rank %d, waiting for site %3d, from %d, to %d... lsite: %3d",
					g.my_rank, glex, source, g.my_rank, lsite);	fflush(stdout);
			}

			//receive data from source

			MPI_Recv( buf, l->num_lattice_site_var, MPI_COMPLEX_PRECISION, source, glex, g.comm_cart, &s);
			//copy received data to output
			vector_PRECISION_copy( phi_out + lsite*l->num_lattice_site_var, buf, 0, l->num_lattice_site_var, l);
			visited_sites[ lsite ] = 1;

			if (g.my_rank == ra) printf(" received!\n");
		    } 

		}
	    }
	}
    }

//    printf("process finished: %d\n", g.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    //copy phi_out to phi
    vector_PRECISION_copy( phi, phi_out, 0, l->num_lattice_site_var * l->num_inner_lattice_sites, l);
  
    int c = 0;
    for (int lsite = 0; lsite < l->num_inner_lattice_sites; lsite++){
	if (visited_sites[lsite] == 0){
	    printf("site %3d not visited on process %d\n", lsite, g.my_rank);
	    c = 1;
	}	
    }
    if (c != 0) exit(0);   
    
   
    
    FREE(phi_out, complex_PRECISION, l->inner_vector_size);
    FREE(visited_sites, int, l->inner_vector_size);
    FREE(buf, complex_PRECISION, l->num_lattice_site_var);
}

#endif
#endif
