#include <stdio.h>
#include "mpi.h"
#include "../MUMPS_5.4.0/include/dmumps_c.h"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

int main(int argc, char ** argv){
	DMUMPS_STRUC_C id;

	int n = 1000;
	int64_t nnz = 1000 + 2*999; // first diag above and below main diag
	int irn[2998];
	int jcn[2998];
	double a[2998];
	double rhs[1000];

	int myid, ierr;
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	
	int i, j;
	irn[0] = 1;
	irn[1] = 2;
	jcn[0] = 1;
	jcn[1] = 1;
	a[0] =  2;
	a[1] = -1;	

		/*
		 2 -1
		-1  2 -1
		   -1  2 -1
			...
				-1  2 -1
				   -1  2
		*/

	rhs[0] = 10;
	rhs[999] = 0;

	irn[2996] =  999;
	irn[2997] = 1000;
	jcn[2996] = 1000;
	jcn[2997] = 1000;
	a[2996] = -1;
	a[2997] =  2;

	for (i = 1; i < 999; i++){
		jcn[3*i -1] = i+1;
		jcn[3*i   ] = i+1;
		jcn[3*i +1] = i+1;
		irn[3*i -1] = i  ;
		irn[3*i   ] = i+1;
		irn[3*i +1] = i+2;
		a[3*i -1] = -1;
		a[3*i   ] =  2;
		a[3*i +1] = -1;
		rhs[i] = 0;
	}


	id.job = JOB_INIT;
	id.par = 1;
	id.sym = 0;
	id.comm_fortran = USE_COMM_WORLD;
	dmumps_c(&id);

	if(myid == 0){
		id.n = n;
		id.nnz = nnz;
		id.irn = irn;
		id.jcn = jcn;
		id.a = a;
		id.rhs = rhs;
	}
#define ICNTL(I) icntl[(I) -1]	//macro according to docu
	id.ICNTL(1) = -1;
	id.ICNTL(2) = -1;
	id.ICNTL(3) = -1;
	id.ICNTL(4) = 0;


	id.job = 6;
	dmumps_c(&id);
	id.job = JOB_END;
	dmumps_c(&id);
	if (myid == 0){
		printf("solution is:\n");
		for (i = 0; i<1000; i++){
			printf("\t\t%8.2f\n", rhs[i]);
		}
	}
	ierr = MPI_Finalize();
	return 0;
}



	
	


















