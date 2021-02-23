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

#ifdef POLYPREC

void local_fgmres_PRECISION_struct_init( local_gmres_PRECISION_struct *p ) {

  p->Z = NULL;
  p->V = NULL;
  p->H = NULL;
  p->x = NULL;
  p->b = NULL;
  p->r = NULL;
  p->w = NULL;
  p->y = NULL;
  p->gamma = NULL;
  p->c = NULL;
  p->s = NULL;
  //p->preconditioner = NULL;
  p->eval_operator = NULL;
}

void local_fgmres_PRECISION_struct_alloc( int m, int n, long int vl, PRECISION tol, const int type, const int prec_kind,
                                          void (*precond)(), void (*eval_op)(), local_gmres_PRECISION_struct *p, level_struct *l ) {


  long int total=0; 
  int i, k=0;
  
  p->restart_length = m;
  p->num_restart = n;

  //p->preconditioner = precond;

  p->eval_operator = eval_op; 
  p->tol = tol;
  p->kind = prec_kind;

#ifdef HAVE_TM1p1
  //vl*=2;
#endif
  
  if(m > 0) {
  total += (m+1)*m; // Hessenberg matrix
  MALLOC( p->H, complex_PRECISION*, m );
  
  total += (5+m)*vl; // x, r, b, w, V
  MALLOC( p->V, complex_PRECISION*, m+1 );

  // no preconditioner
  k = 0;

  total += 4*(m+1); // y, gamma, c, s
  
  p->H[0] = NULL; // allocate connected memory
  MALLOC( p->H[0], complex_PRECISION, total );
  
  p->total_storage = total;
  total = 0;
  
  // ordering: H, y, gamma, c, s, w, V, Z, x, r, b
  // H
  for ( i=1; i<m; i++ )
    p->H[i] = p->H[0] + i*(m+1);
  total += m*(m+1);
  
  // y
  p->y = p->H[0] + total; total += m+1;
  // gamma
  p->gamma = p->H[0] + total; total += m+1;
  // c
  p->c = p->H[0] + total; total += m+1;
  // s
  p->s = p->H[0] + total; total += m+1;
  // w
  p->w = p->H[0] + total; total += vl;
  // V
  for ( i=0; i<m+1; i++ ) {
    p->V[i] = p->H[0] + total; total += vl;
  }
  // Z
  for ( i=0; i<k; i++ ) {
    p->Z[i] = p->H[0] + total; total += vl;
  }

  // x
  p->x = p->H[0] + total; total += vl;
  // r
  p->r = p->H[0] + total; total += vl;
  // b
  p->b = p->H[0] + total; total += vl;
  
  ASSERT( p->total_storage == total );
  }
  
  if ( type == _GLOBAL_FGMRES ) {    
    p->timing = 1;
    p->print = g.vt.evaluation?0:1;
    p->initial_guess_zero = 1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(g.op_PRECISION);
  } else if ( type == _K_CYCLE ) {
    // these settings also work for GMRES as a smoother
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    p->op = &(l->s_PRECISION.op);
  } else if ( type == _COARSE_GMRES ) {
    p->timing = 0;
    p->print = 0;
    p->initial_guess_zero = 1;
    p->layout = -1;
    p->v_start = 0;
    p->v_end = l->inner_vector_size;
    if ( g.odd_even )
      p->op = &(l->oe_op_PRECISION);
    else  
      p->op = &(l->s_PRECISION.op);
  } else {
    ASSERT( type < 3 );
  }
}


void local_fgmres_PRECISION_struct_free( local_gmres_PRECISION_struct *p, level_struct *l ) {

  int k=0;

  //int m = p->restart_length;
  // no preconditioner
  k = 0;

  if(p->restart_length > 0) {
  FREE( p->H[0], complex_PRECISION, p->total_storage );
  FREE( p->H, complex_PRECISION*, p->restart_length );
  FREE( p->V, complex_PRECISION*, p->restart_length+1 );
  
  if ( p->Z != NULL )
    FREE( p->Z, complex_PRECISION*, k );
  }
  
  p->D = NULL;
  p->clover = NULL;
}


void local_process_multi_inner_product_PRECISION( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION psi,
                                                  int start, int end, level_struct *l, struct Thread *threading ) {

  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  int thread_start=start;
  int thread_end=end;

#ifdef _M10TV
  //compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 20);
  for(int c=0; c<count; c++)
    for ( i=thread_start; i<thread_end; )
      FOR20( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
#else
  //compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 2);
  for(int c=0; c<count; c++)
    for ( i=thread_start; i<thread_end; )
      FOR2( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
#endif
}


int local_fgmres_PRECISION( local_gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

  // start and end indices for vector functions, this is always process-specific
  int start;
  int end;

  int j=-1, finish=0, iter=0, il;
  complex_PRECISION gamma0 = 0;

  PRECISION norm_r0=1, gamma_jp1=1;

  if ( l->depth > 0 ) p->timing = 1;

  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  //compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
  start = p->v_start;
  end = p->v_end;

  // towards main loop

  // always assume initial guess set to zero
  vector_PRECISION_copy( p->r, p->b, start, end, l );

  //gamma0 = global_norm_PRECISION( p->r, p->v_start, p->v_end, l, threading ); // gamma_0 = norm(r)
  gamma0 = 0;
  //compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  VECTOR_FOR( int i=start, i<end, gamma0 += NORM_SQUARE_PRECISION(p->r[i]), i++, l );

  p->gamma[0] = gamma0;
  norm_r0 = creal(p->gamma[0]);
  vector_PRECISION_real_scale( p->V[0], p->r, 1/p->gamma[0], start, end, l ); // v_0 = r / gamma_0

  for( il=0; il<p->restart_length && finish==0; il++) {

    j = il; iter++;

    if ( !local_arnoldi_step_PRECISION( p->V, p->Z, p->w, p->H, p->y, j, NULL, p, l, threading ) ) {
      printf0("| -------------- iteration %d, restart due to H(%d,%d) < 0 |\n", iter, j+1, j );
      break;
    }

    if ( cabs( p->H[j][j+1] ) > p->tol/10 ) {
      qr_update_PRECISION( p->H, p->s, p->c, p->gamma, j, l, threading );
      gamma_jp1 = cabs( p->gamma[j+1] );

      if( gamma_jp1/norm_r0 < p->tol || gamma_jp1/norm_r0 > 1E+5 ) { // if satisfied ... stop
        finish = 1;
        if ( gamma_jp1/norm_r0 > 1E+5 ) printf0("Divergence of fgmres_PRECISION, iter = %d, level=%d\n", iter, l->level );
      }
    } else {
      printf0("from fgmres : depth: %d, iter: %d, p->H(%d,%d) = %+lf+%lfi\n", l->depth, iter, j+1, j, CSPLIT( p->H[j][j+1] ) );
      finish = 1;
      break;
    }
  } // end of a single restart

  //compute_solution_PRECISION( p->x, (p->preconditioner&&p->kind==_RIGHT)?p->Z:p->V,
  //                            p->y, p->gamma, p->H, j, (res==_NO_RES)?ol:1, p, l, threading );

  return 0;
}


int local_arnoldi_step_PRECISION( vector_PRECISION *V, vector_PRECISION *Z, vector_PRECISION w,
                                  complex_PRECISION **H, complex_PRECISION* buffer, int j, void (*prec)(),
                                  local_gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

  int i;
  // start and end indices for vector functions depending on thread
  int start, end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  //compute_core_start_end(p->v_start, p->v_end, &start, &end, l, threading);
  start = p->v_start;
  end = p->v_end;

  // no preconditioner
  //apply_operator_PRECISION( w, V[j], p, l, threading ); // w = D*V[j]
  p->eval_operator( w, V[j], p->op, l, threading );

  // orthogonalization
  complex_PRECISION tmp[j+1];
  local_process_multi_inner_product_PRECISION( j+1, tmp, V, w, start, end, l, threading );

  for( i=0; i<=j; i++ )
    H[j][i] = tmp[i];

  for( i=0; i<=j; i++ )
    vector_PRECISION_saxpy( w, w, V[i], -H[j][i], start, end, l );

  // normalization

  PRECISION tmp2 = 0;
  //compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  //PRECISION tmp2 = global_norm_PRECISION( w, p->v_start, p->v_end, l, threading );
  VECTOR_FOR( int i=start, i<end, tmp2 += NORM_SQUARE_PRECISION(p->r[i]), i++, l );

  H[j][j+1] = tmp2;

  // V_j+1 = w / H_j+1,j
  if ( cabs_PRECISION( H[j][j+1] ) > 1e-15 )
    vector_PRECISION_real_scale( V[j+1], w, 1/H[j][j+1], start, end, l );

  return 1;
}


void local_set_ghost_PRECISION( vector_PRECISION phi, const int mu, const int dir,
                                comm_PRECISION_struct *c, const int amount, level_struct *l ) {

  // does not allow sending in both directions at the same time
  if( l->global_splitting[mu] > 1 ) {
    
    int i, j, *table=NULL, mu_dir = 2*mu-MIN(dir,0), offset = c->offset,
        length[2] = {0,0}, comm_start = 0, table_start = 0;
    vector_PRECISION buffer, phi_pt;
    
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir]+c->num_even_boundary_sites[mu_dir]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }
    
#ifdef HAVE_TM1p1
    if ( g.n_flavours == 2 ) {
      length[0] *= 2;
      length[1] *= 2;
      comm_start *= 2;
      offset *= 2;
    }
#endif

    if ( MAX(length[0],length[1]) > c->max_length[mu] ) {
      printf("CAUTION: my_rank: %d, not enough comm buffer\n", g.my_rank ); fflush(0);
      ghost_free_PRECISION( c, l );
      ghost_alloc_PRECISION( MAX(length[0],length[1]), c, l );
    }
    
    buffer = (vector_PRECISION)c->buffer[mu_dir];

    // dir = senddir
    if ( dir == 1 ) {
      memset( buffer, 0.0, length[1]*sizeof(complex_PRECISION) );
    } else if ( dir == -1 ) {
      phi_pt = phi + comm_start;
      memset( phi_pt, 0.0, length[0]*sizeof(complex_PRECISION) );
    } else ASSERT( dir == 1 || dir == -1 );

    // this second part mimics the ghost_wait
    if ( dir == 1 ) {
      int num_boundary_sites = length[0]/offset;

      buffer = (vector_PRECISION)c->buffer[mu_dir];      
      table = c->boundary_table[2*mu+1] + table_start;

      if ( l->depth == 0 ) {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset;
          for ( i=0; i<offset; i++ )
            phi_pt[i] = buffer[i];
          buffer += offset;
        }
      } else {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset;
          for ( i=0; i<offset; i++ )
            phi_pt[i] += buffer[i];
          buffer += offset;
        }
      }
    } else if ( dir == -1 ) {
      // do nothing
    } else ASSERT( dir == 1 || dir == -1 );
  }
}


void coarse_local_n_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                            const int amount, level_struct *l, struct Thread *threading ){

  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      num_4link_var=4*4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      num_link_var=4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;
  config_PRECISION D_pt;

  int core_start;
  int core_end;

  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  //set_boundary_PRECISION( out, 0, l, threading );
  memset( out, 0.0, (l->vector_size-l->inner_vector_size)*sizeof(complex_PRECISION) );

  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  // no cores here, everything is per-process
  core_start = start;
  core_end = start+num_lattice_sites;

  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 0*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }

  // instead of ghost exchanges, set &(op->c) to zero
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      local_set_ghost_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      local_set_ghost_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );    
    }
  }

  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  // no cores here, everything is per-process
  core_start = start;
  core_end = start+num_lattice_sites;

  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }

  END_NO_HYPERTHREADS(threading)
}


void coarse_local_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                          const int amount, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      num_4link_var=4*4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      num_link_var=4*l->num_parent_eig_vect*l->num_parent_eig_vect,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;
  config_PRECISION D_pt;

  int core_start;
  int core_end;

  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  //set_boundary_PRECISION( out, 0, l, threading );
  memset( out, 0.0, (l->vector_size-l->inner_vector_size)*sizeof(complex_PRECISION) );

  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  // no cores here, everything is per-process
  core_start = start;
  core_end = start+num_lattice_sites;

  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 0*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }

  // instead of ghost exchanges, set &(op->c) to zero
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      local_set_ghost_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      local_set_ghost_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );    
    }
  }

  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  // no cores here, everything is per-process
  core_start = start;
  core_end = start+num_lattice_sites;

  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }

  END_NO_HYPERTHREADS(threading)
}


void coarse_local_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in,
                                                    operator_PRECISION_struct *op, level_struct *l,
                                                    struct Thread *threading ) {

  // start and end indices, local i.e. per-process
  int start = op->num_even_sites*l->num_lattice_site_var;
  int end = l->inner_vector_size;

  vector_PRECISION *tmp = op->buffer;
  coarse_diag_ee_PRECISION( out, in, op, l, threading );
  vector_PRECISION_define( tmp[0], 0, start, end, l );

  coarse_local_hopping_term_PRECISION( tmp[0], in, op, _ODD_SITES, l, threading );
  coarse_diag_oo_inv_PRECISION( tmp[1], tmp[0], op, l, threading );
  coarse_local_n_hopping_term_PRECISION( out, tmp[1], op, _EVEN_SITES, l, threading );
}

#endif
