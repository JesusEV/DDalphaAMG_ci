#include <stdio.h>
#include "mpi.h"
#include "/home/leemhuis/installs/MUMPS_5.4.0/include/dmumps_c.h"

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654



#define ICNTL(I) icntl[(I) -1]	//macro according to docu
//#define INFO(I) info[(I) -1]	// does not work, use info(22) instead of info(23)

int main(int argc, char ** argv){
	DMUMPS_STRUC_C id;

	int n = 1000;
	int64_t nnz_loc = 500; // first diag above main diag
	int irn_loc[500];
	int jcn_loc[500];
	double a_loc[500];



	int myid, ierr;
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	

	int i, j;
	

	
	for (i = 0; i < 500; i++){
		irn_loc[i] = myid*500 + i+1;
		jcn_loc[i] = ((myid*500 + i +1)%1000) +1;
		a_loc[i] = 1;

//		rhs[2*i  ] = 1;
//		rhs[2*i+1] = 1;
	}
	

/*
p0	 0 1
  	 0 0 1

		...
		0 0 1
===========================
p1		  0 0 1
		    0 0 1
			...
			  0 0 1
	1 0 .....	... 0 0

		*/

	id.job = JOB_INIT;
	id.par = 1;
	id.sym = 0;
	id.comm_fortran = USE_COMM_WORLD;
	dmumps_c(&id);



	id.ICNTL(5) = 0;
	id.ICNTL(18) = 3;

	id.ICNTL(20) = 11;	//distributed RHS. compare to inctl(20) = 11


	if (myid == 0){
		id.n = n;
	}
	id.nnz_loc = nnz_loc;
	id.irn_loc = irn_loc;
	id.jcn_loc = jcn_loc;
	id.a_loc = a_loc;



//outputs
	id.ICNTL(1) = 6;
	id.ICNTL(2) = -1;
	id.ICNTL(3) = 6;
	id.ICNTL(4) = 2;


//	id.job = 6;	//analyze factorize solve
	id.job = 4;	//analyze factorize
	dmumps_c(&id);





	double rhs_loc[500];
//	double rhs[1000];
	int irhs_loc[500];
	int Nloc_RHS = 500; //no of rows in local rhs
	int LRHS_loc = 500; //leading dimension

	for (i = 0; i < 500; i++){
		rhs_loc[i] = (double)(myid * 500 + i);
		irhs_loc[i] = myid * 500 + i+1;
	}

	if (myid == 0){
		rhs_loc[0] = 77;
		for (i = 0; i < 500; i++){
			printf("%3d\trhs: %5f\n", i+1, rhs_loc[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 1){
		for (i = 0; i < 500; i++){
			printf("%3d\trhs: %5f\n", 501+i, rhs_loc[i]);
		}
	}

	id.nloc_rhs = Nloc_RHS;
	id.rhs_loc = rhs_loc;
	id.irhs_loc = irhs_loc;	
	id.lrhs_loc = LRHS_loc; //leading dimension



	id.job = 9;	//prepare solve with distributed rhs
	

	id.ICNTL(21) = 1;	//non centralized solution
//	if (myid == 0){
//		id.rhs = rhs;
//	}
//	dmumps_c(&id);

	/*
	MPI_Barrier(MPI_COMM_WORLD);	
	printf("myid: %d, info(23): %d\n", myid, id.info[22]);
	printf("\n");
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 0){
		for (i = 0; i < id.info[22]; i++){
			printf("i: %3d, irhs_loc: %d\n", i, irhs_loc[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	exit(0);
	*/



	double SOL_loc[1000];
	int LSOL_loc = id.info[22];
	int ISOL_loc[1000];



	id.sol_loc = SOL_loc;
	id.isol_loc = ISOL_loc;
	id.lsol_loc = LSOL_loc;



	id.job = 3;		// solve
	dmumps_c(&id);
	
	if (myid == 0){
		printf("p0 solution is:\n");
		for (i = 0; i<id.info[22]; i++){
			printf("%4d\t\t, sol:%8.2f, isol: %5d\n", i, SOL_loc[i], ISOL_loc[i]);
		}
	}
/*
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 1){
		printf("\n######################\n");
		for (i = 0; i<1000; i++){
			printf("%4d\t\t%8.2f\n", i, rhs[i]);
		}
	}
*/
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 1){
		printf("p1 solution is:\n");
		for (i = 0; i<id.info[22]; i++){
			printf("%4d\t\t, sol: %8.2f, isol: %5d\n", i, SOL_loc[i], ISOL_loc[i]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);



	id.job = JOB_END;
	dmumps_c(&id);

	printf("process %d finished.\n", myid);
	MPI_Barrier(MPI_COMM_WORLD);



	ierr = MPI_Finalize();
	return 0;
}



	
	


















