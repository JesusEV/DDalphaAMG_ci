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
void blacs_get_(const int*, const int*, int*);

void blacs_pinfo_(int*, int*);
void blacs_gridinit_(int*, const char*, const int*, const int*);
void blacs_gridmap_(int *ictxt, int *usermap, int *ldumap, int *nprow, int *npcol);
void blacs_gridinfo_(const int*, int*, int*, int*, int*);
void descinit_(int*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, int*);
int numroc_(const int*, const int*, const int*, const int*, const int*);

void pgesv_PRECISION(const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, int*, complex_PRECISION*, const int*, const int*, const int*, int* );
void pgetrf_PRECISION(const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, int*, int* );
//		void(const int *, const int *, _Complex float *, const int *, const int *, const int *, int *, int *)
void pgetrs_PRECISION(const char*, const int*, const int*, const complex_PRECISION*, const int*, const int*, const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, int* );
//		 void(const char *, const int *, const int *, const _Complex float *, const int *, const int *, const int *, const int *, _Complex float *, const int *, const int *, const int *, int *)
void pgetri_PRECISION( const int*, complex_PRECISION*, const int*, const int*, int*, int*, complex PRECISION*, const int*, const int*, int*, int* );
void pgemv_PRECISION(char*, int*, int*, PRECISION*, complex_PRECISION*, int*, int*, int*, complex_PRECISION*, int*, int*, int*, int*, PRECISION*, complex_PRECISION*, int*, int*, int*, int*);
void pgemr2d_PRECISION( const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, int*);

void paxpy_PRECISION( const int*, complex_PRECISION*, complex_PRECISION*, const int*, const int*, const int*, const int*, complex_PRECISION*, const int*, const int*, const int*, const int*);
void pnrm2_PRECISION( const int*, PRECISION*, complex_PRECISION*, const int*, const int*, const int*, const int*);
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
      *(g.ds.mumps_Is +k) = i_start + j * site_var + (int)(i/site_var);	// col indices
      *(g.ds.mumps_Js +k) = j_start + j * site_var + (i % site_var); 	// row indices
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
        g.ds.mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] = *(clover_pt + k); //clover is
            // triangular store in DDalphaAMG
        g.ds.mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] = conj_PRECISION(*(clover_pt + k));
      }
      g.ds.mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] = *(clover_pt + k); // diagonal element
    }
    clover_pt += clover_step_size1; // bend pointer to next piece of memory/clover part
    // D
    for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
      for (c = num_eig_vect; c < r; c++, k++){
        g.ds.mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] = *(clover_pt + k);
        g.ds.mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] = conj_PRECISION(*(clover_pt + k));
      }
      g.ds.mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] = *(clover_pt + k);
    }
    clover_pt += clover_step_size1;
    // C
    for (r = num_eig_vect, k = 0; r < 2*num_eig_vect; r++){
      for (c = 0; c < num_eig_vect; c++, k++){
        g.ds.mumps_vals[j*9*SQUARE(site_var) + (r * site_var) + c] = -1.0*(conj_PRECISION(*(clover_pt + k)));
      }
    }

    // no clover_pt correction, clover is hermitian, use same data as for C 
    // B store column-wise / transposed from former storage
    for (r = 0, k = 0; r < num_eig_vect; r++){
      for (c = 0; c < num_eig_vect; c++, k++){
        g.ds.mumps_vals[j*9*SQUARE(site_var) + (c * site_var) + r + num_eig_vect] = *(clover_pt + k);
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
        g.ds.mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] += *(tm_block_pt + k);
        g.ds.mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] -= conj_PRECISION(*(tm_block_pt + k));
      }
      g.ds.mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] += *(tm_block_pt + k);
    }
    tm_block_pt += block_step_size;

    // D
    for (k = 0, r = num_eig_vect; r < 2*num_eig_vect; r++, k++){
      for (c = num_eig_vect; c < r; c++, k++){
        g.ds.mumps_vals[j*9*SQUARE(site_var) + r * site_var + c] += *(tm_block_pt + k);
        g.ds.mumps_vals[j*9*SQUARE(site_var) + c * site_var + r] -= conj_PRECISION(*(tm_block_pt + k));
      }
      g.ds.mumps_vals[j*9*SQUARE(site_var) + r * site_var + r] += *(tm_block_pt + k);
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
        *(g.ds.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir + 1)*num_link_var + k) = 
			-1.0 * *(op->D + 	num_4link_var*op->neighbor_table[index] +
			dir*num_link_var + k); 
	//					find correct block row				start of mu- coupling
        *(g.ds.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) = 
			i_start +	num_site_var * op->neighbor_table[index] + 		k%((int)(num_site_var*0.5));
        //	proc start		block row start						fast changing index
      }

      // C
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var
			+ 1 * (int)SQUARE(num_site_var/2) + k); 
        *(g.ds.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5);
      }

      // B
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var
			+ 2 * (int)SQUARE(num_site_var/2) + k); 
        *(g.ds.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + 0;
      }

      // D
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			-1.0 * *(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var
			+ 3 * (int)SQUARE(num_site_var/2) + k); 
        *(g.ds.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			i_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5);
      }


     // computing and storing column indices Js
      if (comm_nr[dir] > 0 && op->neighbor_table[index+1+dir] >= l->num_inner_lattice_sites){

	li = boundary_table[op->neighbor_table[index + 1 + dir] % comm_nr[dir]]; //will hold the index of neighboring lattice site on neighboring domain

        for (k = 0; k < SQUARE(num_site_var/2); k ++){
//FINDME
          //A
          *(g.ds.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) = 
			neighbors_j_start + num_site_var * li +
			 k/((int)(num_site_var*0.5));
          //C
	  *(g.ds.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * li +
			 k/((int)(num_site_var*0.5));
          //B
          *(g.ds.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * li +
			 k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
          //D
          *(g.ds.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir + 1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) =
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
          *(g.ds.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) =
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5));
          //C 
          *(g.ds.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			 k/((int)(num_site_var*0.5));
          //B
          *(g.ds.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * (op->neighbor_table[index + 1 + dir] %
				l->num_inner_lattice_sites) +
			k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
          //D
          *(g.ds.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) =
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
	  *(g.ds.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + k) = 
			*(buff_d_recv[dir] + buffer_d_pt + k);
	  *(g.ds.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5);
	  *(g.ds.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + k) = 
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5));
        }
	// -C*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(g.ds.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 1 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(g.ds.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5);
	  *(g.ds.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
        }
	// -B*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(g.ds.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 2 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(g.ds.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(g.ds.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5));
        } 
	// D*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(g.ds.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 3 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(g.ds.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(g.ds.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
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
        *(g.ds.mumps_vals + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + k) = 
			-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + k));
        *(g.ds.mumps_Is   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5);
        *(g.ds.mumps_Js   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + k) = 
			(j_start + num_site_var * op->neighbor_table[index] +	 k%(int)(num_site_var*0.5));
      }
      // -C*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds.mumps_vals + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k)); 
        *(g.ds.mumps_Is   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5);
        *(g.ds.mumps_Js   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
      }
      // -B*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds.mumps_vals + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k)); 
        *(g.ds.mumps_Is   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(g.ds.mumps_Js   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5));
      }
      // D*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(g.ds.mumps_vals + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k)); 
        *(g.ds.mumps_Is   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(g.ds.mumps_Js   + (9 * num_link_var)*op->neighbor_table[index + 1 + dir] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			(j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + (int)(num_site_var*0.5));
      }
    }	//loop over nodes  
    //printf0("mu- on node %d, dir %d done!\n", node, dir);
  }	//loop over directions

 
 //printf0("hopping in mu - done!\n");

#ifdef COARSE_SCALAP	//use Scalapack on coarsest level to solve.
  coarse_scalap_setup_PRECISION( l, threading);
 //printf0("generate scalap matrix done\n");
#else
  // increase global indices by 1 to match fortran indexing.
  // spmv doesn't work then anymore
  int nnz_loc = SQUARE(site_var) * nr_nodes *9;
  for (i = 0; i < nnz_loc; i++){	//increase indices by one to match fortran indexing in MUMPS
    *(g.ds.mumps_Js + i ) = *(g.ds.mumps_Js + i ) +1;
    *(g.ds.mumps_Is + i ) = *(g.ds.mumps_Is + i ) +1;
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

void direct_solves_set_reset_PRECISION( level_struct *l, struct Thread *threading ){
  g.coarsest_time = 0;
    if (g.on_solve) {
      for (int i = 1; i<g.num_levels; i++) l = l->next_level;
      if (!l->idle){ 

	  printf0("call to setup from top_level.c\n");
	  mumps_setup_PRECISION(l, threading);        //setup vals, Is, Js

          double t0,t1;
          START_MASTER(threading)
          t0 = MPI_Wtime();

#if defined(MUMPS_ADDS)
          printf0("starting analyze from top_level.c\n");

	  //    g.mumps_id.job = 4; //analyze and factorize
          g.mumps_id.job = 1; //analyze

          cmumps_c(&(g.mumps_id));

          printf0("analyze done, starting factorize singlethreaded from top_level.c\n");
          g.mumps_id.job = 2; //factorize
          cmumps_c(&(g.mumps_id)); //only factorize when on solve
#elif defined(COARSE_SCALAP)
	  coarse_scalap_factorize_PRECISION( l, g.ds.dense_vals2d,
		  g.ds.desc_dense_vals2d, g.ds.ipiv, threading );//only factorize when on solve
#endif

	  t1 = MPI_Wtime();
#if defined(MUMPS_ADDS)
	  printf0("MUMPS analyze and factorize time (seconds) : %f \t from top_level.c\n",t1-t0);
#elif defined(COARSE_SCALAP)
	  printf0("Invert using scalapack time (seconds) : %f \t from top_level.c\n",t1-t0);
#endif
	  g.coarsest_fact_time += t1-t0;
	  END_MASTER(threading)
	  SYNC_CORES(threading)

	  //set direct solve as precond during solve phase
#if defined(MUMPS_ADDS)        
          l->p_PRECISION.preconditioner = mumps_solve_PRECISION;
	  l->p_PRECISION.eval_operator = coarse_apply_oddeven_operator_PRECISION;
#endif
	  //use entire coarsest system, not only odd-even Schur-Decomp when using direct solves
	  if (g.odd_even) l->p_PRECISION.v_end *=2;
        }
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
        *(g.ds.mumps_irhs_loc + i) = lrank * rhs_len + i+1;		//+1 due to fortran indexing
      }

      // copying eta to local_rhs
      vector_PRECISION_copy(g.ds.mumps_rhs_loc, eta, px->v_start, px->v_end, lx );

      // centralized solution
      if (g.my_rank == 0){
	// FIXME : do some sort of casting here, to avoid warnings at compile-time
	g.mumps_id.rhs = (mumps_complex *)g.ds.mumps_SOL;
      }

      // solving
      g.mumps_id.job = 3; // solve
      cmumps_c(&(g.mumps_id));

      // distributing the solution to all processes. Must be stored in px->x
      int send_count = (lx->p_PRECISION.v_end-lx->p_PRECISION.v_start);
      MPI_Scatter(g.ds.mumps_SOL, send_count, MPI_COMPLEX_PRECISION, phi, send_count, MPI_COMPLEX_PRECISION, 0, lx->gs_PRECISION.level_comm); // lx->gs_PRECISION.level_comm

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
void coarse_scalap_solve_PRECISION(vector_PRECISION phi, vector_PRECISION Dphi,
                           vector_PRECISION eta, int res, level_struct *l,
                           struct Thread *threading){
    if (!l->idle){
//	printf0("solving scalap...\n");
        START_MASTER(threading)
        g.coarsest_solve_time -= MPI_Wtime();
	END_MASTER(threading)
        SYNC_CORES(threading);
	
	int ione = 1, info = 0;
	char trans = 'N';

	int N = l->num_processes * l->num_inner_lattice_sites * l->num_lattice_site_var;

	PRECISION alpha = 1.0, beta = 0.0;

	
	//rhs2d = eta
	scalap_1d_2d_vec_PRECISION( eta, l, threading);
	/*
	rhs2d = (LU)^-1 rhs2d
	pgetrs_PRECISION( &trans, &N, &ione, g.ds.dense_vals2d, &ione, &ione,
		g.ds.desc_dense_vals2d, g.ds.ipiv,
		g.ds.rhs2d, &ione, &ione, g.ds.desc_rhs2d, &info );
		*/
	pgemv_PRECISION( &trans, &N, &N, &alpha, g.ds.dense_vals2d, &ione, &ione,
		g.ds.desc_dense_vals2d, g.ds.rhs2d, &ione, &ione,
		g.ds.desc_rhs2d, &ione, &beta, g.ds.rhs2d, &ione, &ione,
		g.ds.desc_rhs2d, &ione); 
	if (info != 0 ) error0("Error during pgetrs_(), info = %d\n", info);
	
	//phi = rhs2d (which is sol)
	scalap_2d_1d_vec_PRECISION( phi, l, threading);

        g.coarsest_solve_number ++;
        g.coarsest_solve_time += MPI_Wtime();
        printf0("scalap time  = %f, scalap solves:  %d\n", g.coarsest_solve_time, g.coarsest_solve_number);
    }
}


void coarse_scalap_factorize_PRECISION( level_struct *l, vector_PRECISION A, int* descA, int* ipiv, struct Thread *threading){
   
    printf0("scalapack factorization ongoing...\n");
    fflush(stdout);
   //distribute Matrix to 2d cyclic pattern
    scalap_1d_2d_A_PRECISION( l, threading);
   
    int info = 0;
    int ione = 1; //starting indices (global)
    int N = l->num_inner_lattice_sites * l->num_lattice_site_var * l->num_processes; 
       //		    ( M,    N,	    A,		IA, JA, DESCA, IPIV, INFO)
    if ( g.ds.myrow2d >= 0){
	pgetrf_PRECISION( &N, &N, g.ds.dense_vals2d, &ione, &ione, g.ds.desc_dense_vals2d, ipiv, &info );    
        if (info != 0 ) error0("Error during pgetrf_(), info = %d\n", info);
    
	complex PRECISION work_query;
	int iwork_query;

	int lwork = -1, liwork = -1;
	pgetri_PRECISION(
	    &N,
	    g.ds.dense_vals2d,              /* LU-Faktorisierte Matrix */
	    &ione, &ione,           /* IA, JA */
	    g.ds.desc_dense_vals2d,
	    ipiv,
	    &work_query,
	    &lwork,
	    &iwork_query,
	    &liwork,
	    &info
	);
        if (info != 0 ) error0("Error during pgetri_query(), info = %d\n", info);
	
	complex_PRECISION *work = NULL;
	int *iwork = NULL;

	lwork = (int) creal(work_query);
	liwork = iwork_query;
	work = malloc(lwork * sizeof(complex_PRECISION));
	iwork = malloc (liwork * sizeof(int));
    	pgetri_PRECISION(
	    &N,
	    g.ds.dense_vals2d,
	    &ione, &ione,
	    g.ds.desc_dense_vals2d,
	    ipiv,
	    work,
	    &lwork,
	    iwork,
	    &liwork,
	    &info
	);
        if (info != 0 ) error0("Error during pgetri_(), info = %d\n", info);
    
	free(work);
	free(iwork);
    
    
    }

}


void coarse_scalap_setup_PRECISION(level_struct *l, struct Thread *threading){

    int r, c; //for row and column index for a given matrix element
    memset(g.ds.dense_vals, 0, l->num_inner_lattice_sites * l->num_lattice_site_var *
	    l->num_inner_lattice_sites * l->num_lattice_site_var * l->num_processes *
	    sizeof(complex_PRECISION));


    int global_comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &global_comm_size);

    int *glob_ranks = NULL;
    int *loc_ranks = NULL; //will contain global ranks and corresponding local ranks

    MALLOC( glob_ranks, int, global_comm_size);
    MALLOC( loc_ranks, int, global_comm_size);
    for (int i = 0; i<global_comm_size; i++) { glob_ranks[i]= i; loc_ranks[i] = -1;}
    
    MPI_Group world_group;
    MPI_Group local_group;

    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Comm_group(l->gs_PRECISION.level_comm, &local_group);

    MPI_Group_translate_ranks(world_group, global_comm_size, glob_ranks, local_group, loc_ranks);

  //translates from global to local ranks


    for (int llsite = 0; llsite < l->num_inner_lattice_sites; llsite++ ){ //loop over lattice sites
    //each lattice site contains 9 blocks of each SQUARE(site_var) elements, 1 for self coupling, +
    //2x4 hopping terms (2 in each dimension)
	for (int d = 0; d < 9; d++){ //loop over each block as mentioned above
	    for (int i = 0; i < l->num_lattice_site_var; i++)
		for (int j = 0; j < l->num_lattice_site_var; j++){ //local index to copy elementwise within a block
		r = g.ds.mumps_Is[llsite*9*SQUARE(l->num_lattice_site_var) + d * SQUARE(l->num_lattice_site_var) + 
		    i * l->num_lattice_site_var + j] - 
			   loc_ranks[g.my_rank] * l->num_inner_lattice_sites * l->num_lattice_site_var;
		c = g.ds.mumps_Js[llsite*9*SQUARE(l->num_lattice_site_var) + d * SQUARE(l->num_lattice_site_var) + 
		    i * l->num_lattice_site_var + j];
		g.ds.dense_vals[c * l->num_inner_lattice_sites * l->num_lattice_site_var +
		    r] += g.ds.mumps_vals[llsite*9*SQUARE(l->num_lattice_site_var) + d *
		    SQUARE(l->num_lattice_site_var) + i * l->num_lattice_site_var + j];
		}    
	}
    }
   
    free(glob_ranks);
    free(loc_ranks);

}


void coarse_scalap_init_PRECISION(level_struct *l, struct Thread *threading){


    //1. Initialize 1D - (cyclic) Blacs-Grid for reading in/out data
    int izero = 0;
    int ione = 1;
    int nprow = l->num_processes;
    int npcol = 1;
    int info = 0;
    int mycol;
    int N = l->num_processes * l->num_inner_lattice_sites * l->num_lattice_site_var;

    int ictxt1d = MPI_Comm_c2f(l->gs_PRECISION.level_comm);
       
    int usermap[4] = {0, 1, 2, 3};
    blacs_gridmap_( &ictxt1d, usermap, &nprow, &nprow, &npcol ); 
    blacs_gridinfo_(&ictxt1d, &nprow, &npcol, &g.ds.myrow, &mycol);	//set process coordinates of blacs-grid
    
    g.ds.blacs_ctxt1d = ictxt1d; //blacs_gridmap_ changes ictxt, therefore save new context


    //setting the descriptors:
    int bs = l->num_lattice_site_var * l->num_inner_lattice_sites;
    int nrhs = 1;
	   
    int numr = numroc_( &N, &bs, &g.ds.myrow, &izero, &nprow ); // number of rows stored in each process
    int lldA = numr > 1? numr : 1;	//leading dimension in A (remember, matrix elements are	stored in a column major order)

    //descinit ( DESC,			    M, N,    MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO )
    descinit_( g.ds.desc_dense_vals, &N, &N, &bs, &bs, &izero, &izero, &ictxt1d, &lldA, &info);
    if (info != 0) error0("Error in descinit for DescA, info = %d\n", info);

    //descinit ( DESC,		    M, N,	MB, NB,	    IRSRC, ICSRC, ICTXT, LLD, INFO )
    descinit_( g.ds.desc_rhs, &N, &nrhs, &bs, &ione, &izero, &izero, &ictxt1d, &lldA, &info);
    if (info != 0) error0("Error in descinit for DescB, info = %d\n", info);

    //2. Initialize 2D - cyclic Blacs-Grid for performing fast calculation
    int nprow2d = g.prow2d;
    int npcol2d = g.pcol2d;
    int mycol2d;


    //TODO: Remove this one the second comm is MPI_COMM_WORLD:	<-- remove this line!
//    MPI_Comm level_comm_dup;
//    MPI_Comm_dup(l->gs_PRECISION.level_comm, &level_comm_dup);

//    int ictxt2d = MPI_Comm_c2f(level_comm_dup);
    printf("r: %d, CHECKPOINT1\n", g.my_rank); fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    int ictxt2d = MPI_Comm_c2f(MPI_COMM_WORLD);
    
    blacs_gridmap_( &ictxt2d, usermap, &nprow2d, &nprow2d, &npcol2d ); 
    blacs_gridinfo_( &ictxt2d, &nprow2d, &npcol2d, &g.ds.myrow2d, &mycol2d );

    g.ds.blacs_ctxt2d = ictxt2d; //blacs_gridmap_ changes ictxt, therefore save new context

    int bs2d = g.bs2d;
    int numr2d = numroc_( &N, &bs2d, &g.ds.myrow2d, &izero, &nprow2d );
    int lldA2d = numr2d > 1? numr2d : 1;


    if (g.ds.myrow2d >= 0) { 
	 //descinit ( DESC,			    M, N,    MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO )
	descinit_( g.ds.desc_dense_vals2d, &N, &N, &bs2d, &bs2d, &izero, &izero, &ictxt2d, &lldA2d, &info);
	if (info != 0) error0("Error in descinit for DescA2D, info = %d\n", info);

	//descinit ( DESC,		    M, N,	MB, NB,	    IRSRC, ICSRC, ICTXT, LLD, INFO )
	descinit_( g.ds.desc_rhs2d, &N, &nrhs, &bs2d, &ione, &izero, &izero, &ictxt2d, &lldA2d, &info);
	if (info != 0) error0("Error in descinit for DescB2D, info = %d\n", info);
    } else {
	g.ds.desc_dense_vals2d[0] = 1;
	g.ds.desc_dense_vals2d[1] = ictxt2d;
	g.ds.desc_dense_vals2d[2] = 0;
	g.ds.desc_dense_vals2d[3] = 0;
	g.ds.desc_dense_vals2d[4] = 0;
	g.ds.desc_dense_vals2d[5] = 0;
	g.ds.desc_dense_vals2d[6] = 0;
	g.ds.desc_dense_vals2d[7] = 0;
	g.ds.desc_dense_vals2d[8] = 0;

	g.ds.desc_rhs2d[0] = 1;
	g.ds.desc_rhs2d[1] = ictxt2d;
	g.ds.desc_rhs2d[2] = 0;
	g.ds.desc_rhs2d[3] = 0;
	g.ds.desc_rhs2d[4] = 0;
	g.ds.desc_rhs2d[5] = 0;
	g.ds.desc_rhs2d[6] = 0;
	g.ds.desc_rhs2d[7] = 0;
	g.ds.desc_rhs2d[8] = 0;
    }
}

void scalap_1d_2d_A_PRECISION(level_struct *l, struct Thread *threading){
    int N = l->num_processes * l->num_inner_lattice_sites * l->num_lattice_site_var;
    int ione = 1; //starting indices
   
    //pgemr2d_PRECISION( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
    if (g.ds.myrow >= 0 || g.ds.myrow2d >= 0){
        pgemr2d_PRECISION( &N, &N, g.ds.dense_vals, &ione, &ione, g.ds.desc_dense_vals,
	    g.ds.dense_vals2d, &ione, &ione, g.ds.desc_dense_vals2d,
	    &(g.ds.blacs_ctxt2d));
    }
}

void scalap_2d_1d_A_PRECISION(level_struct *l, struct Thread *threading){
    /*
       probably never used!
       */
    int N = l->num_processes * l->num_inner_lattice_sites * l->num_lattice_site_var;
    int ione = 1; //starting indices
    int ictxt = g.ds.blacs_ctxt2d;

    //pgemr2d_PRECISION( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
    if (g.ds.myrow >= 0 || g.ds.myrow2d >= 0){
	pgemr2d_PRECISION( &N, &N, g.ds.dense_vals2d, &ione, &ione,
	    g.ds.desc_dense_vals2d,
	    g.ds.dense_vals, &ione, &ione, g.ds.desc_dense_vals, &ictxt);
    }
}

void scalap_1d_2d_vec_PRECISION( vector_PRECISION vec, level_struct *l, struct Thread *threading){
    int N = l->num_processes * l->num_inner_lattice_sites * l->num_lattice_site_var;
    int ione = 1; //starting indices
    //(global) The row and column indices in the array A indicating the first row and the first column, respectively, of the submatrix of A) to copy. 
    int ictxt = g.ds.blacs_ctxt2d;

    //pgemr2d_PRECISION( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
    if (g.ds.myrow >= 0 || g.ds.myrow2d >= 0){
	pgemr2d_PRECISION( &N, &ione, vec, &ione, &ione, g.ds.desc_rhs, g.ds.rhs2d, &ione, &ione, g.ds.desc_rhs2d, &ictxt);
    }
}

void scalap_2d_1d_vec_PRECISION( vector_PRECISION vec, level_struct *l, struct Thread *threading){
    int N = l->num_processes * l->num_inner_lattice_sites * l->num_lattice_site_var;
    int ione = 1; //starting indices
    //(global) The row and column indices in the array A indicating the first row and the first column, respectively, of the submatrix of A) to copy. 
    int ictxt = g.ds.blacs_ctxt2d;

    //pgemr2d_PRECISION( m, n, a, ia, ja, desca, b, ib, jb, descb, ictxt);
    if (g.ds.myrow >= 0 || g.ds.myrow2d >= 0){
	pgemr2d_PRECISION( &N, &ione, g.ds.rhs2d, &ione, &ione, g.ds.desc_rhs2d, vec, &ione, &ione, g.ds.desc_rhs, &ictxt);
    }
}
#endif
#endif

