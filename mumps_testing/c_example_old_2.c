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
	int irn_loc[1499];	//1499 * 2 = 2998
	int jcn_loc[1499];
	double a_loc[1499];
	double rhs_loc[500];

	int myid, ierr;
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	

	int i, j;

	if (myid == 0){
		irn_loc[0] = 1;
		irn_loc[1] = 2;
		jcn_loc[0] = 1;
		jcn_loc[1] = 1;
		a_loc[0] =  2;
		a_loc[1] = -1;	

		rhs_loc[0] = 10;

		/*
		 4 -1
		-1  4 -1
		   -1  4 -1
			...
p0			-1  4 -1
=========================================================
p1			   -1  4 -1
				... 
				-1  4 -1
				   -1  4
		*/


		for (i = 1; i < 499; i++){
			irn_loc[3*i -1] = i;
			irn_loc[3*i   ] = i +1;
			irn_loc[3*i +1] = i +2;

			jcn_loc[3*i -1] = i +1;
			jcn_loc[3*i   ] = i +1;
			jcn_loc[3*i +1] = i +1;

			a_loc[3*i -1] = -1;
			a_loc[3*i   ] =  2;
			a_loc[3*i +1] = -1;

			rhs_loc[i] = 0;
		}

		irn_loc[1496] = 499;
		irn_loc[1497] = 500;
		irn_loc[1498] = 500;
					//p1 holds 0 - 1498
					//p2 holds 1499 - 2998
		jcn_loc[1496] = 500;
		jcn_loc[1497] = 500;
		jcn_loc[1498] = 501;
		a_loc[1496] = -1;
		a_loc[1497] =  2;
		a_loc[1498] = -1;
		rhs_loc[499] = 0;
	}
	else if (myid == 1){
		irn_loc[0] = 501;
		irn_loc[1] = 501;
		irn_loc[2] = 502;
		jcn_loc[0] = 500;
		jcn_loc[1] = 501;
		jcn_loc[2] = 501;
		a_loc[0] = -1;
		a_loc[1] =  2;
		a_loc[2] = -1;
		rhs_loc[0] = 0;
		
		for (i = 1; i < 499; i++){
			jcn_loc[3*i   ] = 500 + i;
			jcn_loc[3*i +1] = 500 + i;
			jcn_loc[3*i +2] = 500 + i;
			irn_loc[3*i   ] = 500 + i    ;
			irn_loc[3*i +1] = 500 + i + 1;
			irn_loc[3*i +2] = 500 + i + 2;
	
			a_loc[3*i   ] = -1;
			a_loc[3*i +1] =  2;
			a_loc[3*i +2] = -1;
			rhs_loc[i] = 0;
		}
		irn_loc[1497] = 999;
		irn_loc[1498] = 1000;

		jcn_loc[1497] = 1000;
		jcn_loc[1498] = 1000;
	
		a_loc[1497] = -1;
		a_loc[1498] =  2;

		rhs_loc[499] = 1;
	}


	id.job = JOB_INIT;
	id.par = 1;
	id.sym = 0;
	id.comm_fortran = USE_COMM_WORLD;
	dmumps_c(&id);


#define ICNTL(I) icntl[(I) -1]	//macro according to docu

	id.ICNTL(5) = 0;
	id.ICNTL(18) = 3;

	if (myid == 0){
		id.n = 1000;
	}
	id.nnz_loc = 1499;
	id.irn_loc = irn_loc;
	id.jcn_loc = jcn_loc;
	id.a_loc = a_loc;
	id.rhs_loc = rhs_loc;





//outputs
	id.ICNTL(1) = 6;
	id.ICNTL(2) = -1;
	id.ICNTL(3) = 6;
	id.ICNTL(4) = 2;


	id.job = 6;
	dmumps_c(&id);
	printf("ICNTL(18): %3d\n\n\n", id.ICNTL(18));

	id.job = JOB_END;
	dmumps_c(&id);
	if (myid == 0){
		printf("solution is:\n");
		for (i = 0; i<500; i++){
			printf("%4d\t\t%8.2f\n", i, rhs_loc[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 1){
		for (i = 0; i<500; i++){
			printf("%4d\t\t%8.2f\n", 500+i, rhs_loc[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	printf("process %d finished.\n", myid);
	MPI_Barrier(MPI_COMM_WORLD);

	ierr = MPI_Finalize();
	return 0;
}



	
	


















