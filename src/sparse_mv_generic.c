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
    if( proc_neighbors[i]==g.my_rank ){
      ins[i] = in;
    }
    else{
      ins[i] = NULL;
      MALLOC( ins[i],complex_PRECISION,vl );
      // hardcoded the <5> as a tag, as MPI_ANY_TAG was (weirdly) throwing an error
      //printf("proc=%d, neigh=%d\n", g.my_rank, proc_neighbors[i]);
      MPI_Isend( in, vl, MPI_COMPLEX_PRECISION, proc_neighbors[i], 5, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm, &(reqs[i]) );
      MPI_Irecv( ins[i], vl, MPI_COMPLEX_PRECISION, proc_neighbors[i], 5, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm, &(reqr[i]) );
    }
  }

  for( i=0;i<8;i++ ){
    if( proc_neighbors[i]!=g.my_rank ) {
      MPI_Wait( &(reqs[i]),&stat );
      MPI_Wait( &(reqr[i]),&stat );
    }
  }

  vector_PRECISION_define( out,0,p->v_start,p->v_end,l );

  //printf("%d\n", vl);
  //printf("%d\n", n);
  //printf("%d\n", n/48/g.num_processes);




  /*
  if( g.my_rank==0 ){
    int offst = 2;
    vector_PRECISION Ax = A + (48*48)*offst;
    int* Isx = Is + (48*48)*offst;
    int* Jsx = Js + (48*48)*offst;
  
    for( i=0;i<48*48;i++ ){
      printf("%d,%d\t", Isx[i], Jsx[i]);
      if( i%48==0 ){ printf("\n"); }
    }

    for( i=0;i<48*48;i++ ){
      printf("%f+i%ft", CSPLIT(Ax[i]));
      if( i%48==0 ){ printf("\n"); }
    }
  }
  */




  //MPI_Barrier( MPI_COMM_WORLD );
  //exit(0);


  //int* bfr = Is;
  //Is = Js;
  //Js = bfr;
//  printf("%d", vl);
//  exit(0);


  for( w=0;w<n;w++ ){
    i = Is[w]%vl;
    jx = Js[w];
    //if(g.my_rank==1 && Is[w]!=0) printf("%d\n", l->num_inner_lattice_sites * g.num_processes * l->num_lattice_site_var);
    //if(g.my_rank==1 && Is[w]!=0) printf("%d\n", Is[w]);
    //if(g.my_rank==1 && Is[w]!=0) printf("%d\n\n", Js[w]);
    //continue;
    px = jx / vl;
//    px = jx/(vl*l->num_lattice_site_var);
    if( px==g.my_rank ) inx = in;
    else{
      for( wx=0;wx<8;wx++ ){
        if( px==proc_neighbors[wx] ) inx = ins[wx];
      }
    }

    j = jx%vl;
    //printf("i=%d,j=%d, lth=%d, vl=%d, vl/48=%d, l->num_lattice_site_var=%d, n=%d\n", i,j,vl*l->num_lattice_site_var, vl, vl/48, l->num_lattice_site_var, n);
    //printf("j=%d,jx=%d\n", j,jx);
    //printf("in=%p,inx=%p\n", in,inx);
    rval = inx[j];
    out[i] += A[w]*rval;
  }


  //exit(0);


  END_MASTER(threading)

}
