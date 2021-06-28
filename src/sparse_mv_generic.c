#include "main.h"

/*
input parameters:
  -- out       : the resulting vector in out=A*in
  -- in        : the input vector in out=A*in
  -- A         : the matrix values
  -- Is        : the row indices of the matrix A
  -- Js        : the column indices of the matrix A
  -- n         : the number of nonzero elements in the sparse matrix A
  -- p         : information typically used in FGMRES within DDalphaAMG
  -- l         : level information within DDalphaAMG
  -- threading : threading information within DDalphaAMG
*/


void spmv_PRECISION( vector_PRECISION out, vector_PRECISION in, vector_PRECISION A, int *Is, int *Js,
                     int n, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

  START_MASTER(threading)

  int i, vl = p->v_end-p->v_start, proc_neighbors[8], jx, j, px, w, wx;
  complex_PRECISION rval;
  complex_PRECISION *inx;
  vector_PRECISION ins[9];
  MPI_Request reqs[8],reqr[8];
  MPI_Status stat;

  for( i=0;i<8;i++ ){
    proc_neighbors[i] = l->neighbor_rank[i];
  }

  ins[8] = in;

  for( i=0;i<8;i++ ){
    if( proc_neighbors[i]==g.my_rank ) ins[i] = in;
    else{
      ins[i] = NULL;
      MALLOC( ins[i],complex_PRECISION,vl );
      MPI_Isend( in, vl, MPI_COMPLEX_PRECISION, proc_neighbors[i], MPI_ANY_TAG, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm, &reqs[i] );
      MPI_Irecv( ins[i], vl, MPI_COMPLEX_PRECISION, proc_neighbors[i], MPI_ANY_TAG, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm, &reqr[i] );
    }
  }

  for( i=0;i<8;i++ ){
    if( proc_neighbors[i]!=g.my_rank ) {
      MPI_Wait( &reqs[i],&stat );
      MPI_Wait( &reqr[i],&stat );
    }
  }

  vector_PRECISION_define( out,0,p->v_start,p->v_end,l );

  for( w=0;w<n;w++ ){
    i = Is[w];
    jx = Js[w];
    px = jx/vl;
    if( px==g.my_rank ) inx = in;
    else{
      for( wx=0;wx<8;wx++ ){
        if( px==proc_neighbors[wx] ) inx = ins[wx];
      }
    }
  }

  j = jx%vl;
  rval = inx[j];
  out[i] += A[w]*rval;

  END_MASTER(threading)
}
