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

#ifdef MUMPS_ADDS

#include "mumps_PRECISION.h"


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
#endif

  printf0("twisted mass part done!\n");

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
      
	    
      printf0("starting node: %d, core_end: %d\n", node, core_end);



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
          //A
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + k) = 
			neighbors_j_start + num_site_var * boundary_table[bt_index] + k/((int)(num_site_var*0.5));
          //C
	  *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * boundary_table[bt_index] + k/((int)(num_site_var*0.5));
          //B
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * boundary_table[bt_index] + k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);
          //D
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir + 1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			neighbors_j_start + num_site_var * boundary_table[bt_index] + k/((int)(num_site_var*0.5)) + (int)(num_site_var*0.5);

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
          		j_start + num_site_var * op->neighbor_table[index +1 + dir] + k/((int)(num_site_var*0.5));
          //C 
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir +1)*num_link_var + 1 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * op->neighbor_table[index +1 + dir] + k/(int)(num_site_var*0.5) + 0;
          //B
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 2 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * op->neighbor_table[index +1 + dir] + k/(int)(num_site_var*0.5) + (int)(num_site_var*0.5);
          //D
          *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + (2*dir+1)*num_link_var + 3 * (int)SQUARE(num_site_var/2) + k) = 
			j_start + num_site_var * op->neighbor_table[index +1 + dir] + k/(int)(num_site_var*0.5) + (int)(num_site_var*0.5);
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
			neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5);
        }
	// -C*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_PRECISION.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 1 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(l->p_PRECISION.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5);
	  *(l->p_PRECISION.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) =
			neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5) + num_site_var*0.5;
        }
	// -B*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_PRECISION.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 2 * SQUARE((int)(num_site_var*0.5)) + k); 
	  *(l->p_PRECISION.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(l->p_PRECISION.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5);
        } 
	// D*
	for (k = 0; k < SQUARE(num_site_var / 2); k++ ){
	  *(l->p_PRECISION.mumps_vals + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) =
			*(buff_d_recv[dir] + buffer_d_pt + 3 * SQUARE((int)(num_site_var*0.5)) + k);
	  *(l->p_PRECISION.mumps_Is + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
	  *(l->p_PRECISION.mumps_Js + 9 * num_link_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			neighbors_j_start + num_site_var * buff_i_recv[dir][2 * buffer_i_pt + 1] + k%(int)(num_site_var*0.5) + num_site_var*0.5;
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
			j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5);
      }
      // -C*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5);
        *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 1 * SQUARE((int)(num_site_var*0.5)) + k) = 
			j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + num_site_var*0.5;
      }
      // -B*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 2 * SQUARE((int)(num_site_var*0.5)) + k) = 
			j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5);
      }
      // D*
      for (k = 0; k < SQUARE(num_site_var/2); k ++){
        *(l->p_PRECISION.mumps_vals + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k ) = 
			-1.0 * conj_PRECISION(*(op->D + num_4link_var*op->neighbor_table[index] + dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k));
        *(l->p_PRECISION.mumps_Is + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			i_start + num_site_var * op->neighbor_table[index + 1 + dir] + k/(int)(num_site_var*0.5) + num_site_var*0.5;
        *(l->p_PRECISION.mumps_Js + (9 * num_link_var)*op->neighbor_table[index] + num_link_var + 2*dir*num_link_var + 3 * SQUARE((int)(num_site_var*0.5)) + k) = 
			j_start + num_site_var * op->neighbor_table[index] + k%(int)(num_site_var*0.5) + num_site_var*0.5;
      }
    }	//loop over nodes  
    printf0("mu- on node %d, dir %d done!\n", node, dir);
  }	//loop over directions


  printf0("hopping in mu - done!\n");

   

  // increase global indices by 1 to match fortran indexing.
  // spmv doesn't work then anymore
  int nnz_loc = SQUARE(site_var) * nr_nodes *9;
  for (i = 0; i < nnz_loc; i++){	//increase indices by one to match fortran indexing
    *(l->p_PRECISION.mumps_Js + i ) = *(l->p_PRECISION.mumps_Js + i ) +1;
    *(l->p_PRECISION.mumps_Is + i ) = *(l->p_PRECISION.mumps_Is + i ) +1;
  }

  printf0("freeing ranks\n");
  //TODO: fix threading in here!
  //START_MASTER(threading)
  FREE( glob_ranks, int, global_comm_size);
  FREE( loc_ranks, int, global_comm_size);
//  END_MASTER(threading)
//  SYNC_CORES(threading)



  // timing the setup
//  START_MASTER(threading)
  t1 = MPI_Wtime();
  printf0("MUMPS pre-setup time (seconds) : %f\n",t1-t0);
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  

}



void mumps_solve_PRECISION( vector_PRECISION phi, vector_PRECISION Dphi, vector_PRECISION eta,
                            int res, level_struct *lx, struct Thread *threading )
{
    if (!lx->idle){
      START_MASTER(threading)
      g.mumps_solve_time -= MPI_Wtime();
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
      g.mumps_solve_number ++;
      g.mumps_solve_time += MPI_Wtime();
      printf0("mumps time  = %f, mumps solves:  %d\n", g.mumps_solve_time, g.mumps_solve_number);
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
