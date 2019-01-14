#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include <pthread.h>
#include <mpi.h>


double** allocMat(int size){
  	// Allocation of Matrix
    double **mat = (double **)malloc(size * sizeof(double *)); 
    for (int i=0; i<size; i++)
    	mat[i] = (double *)malloc(size * sizeof(double)); 
   return mat; 
}

void initMat(double** mat, int size, double init_val){
	int i, j; 
    
    for (i = 0; i <  size; i++)
    {
		for (j = 0; j < size; j++)
		{
			if(j+1==size || i+1 == size || j== 0 || i == 0)		mat[i][j] = 0;
			else	mat[i][j] = init_val;
		}
	}
   // return mat; 
}

void initCircle(double** mat, int param_grid_size, int diameter_m, double H)
{
	// int grid_size = param_grid_size;
	double center = float(param_grid_size-1)/2;
	int z,s;
	for (z=1; z<param_grid_size-1; ++z)
	{
		for (s=1; s<param_grid_size-1; ++s)
		{
			if (pow(float(z-center),2) + pow(float(s-center),2) <= pow(float(diameter_m)/2,2))	mat[z][s] = H;
		}
	}
}

void copyMatrix(double** mat, double** mat_backup , int size)
{	
	for(int z=0; z<size; ++z)
	{
		for(int s=0; s<size; ++s)
			mat_backup[z][s] = mat[z][s];
	}
}	



pthread_t* initThreads(int numThreads)
{
	pthread_t* threads =  (pthread_t* )malloc(numThreads * sizeof(pthread_t));
	return(threads);
}

//struct für die Parameterübergabe an Threads
struct args {
		int zeile;
		int spalte;
		double** mat;
		double** mat_backup;
		int size;
		double* G;
		int numThreads;
		~args(){}
	};
 
 void *compute(void *input){
	int z = ((args*)input)->zeile;
	int size = ((args*)input)->size;
	for(int j = 1; j < size-1 ; ++j){
		for(int k = 0; k < 3; ++k){
			((args*)input)->mat[z][j] += ((args*)input)->G[k]*(((args*)input)->mat_backup[z][j+k-1] + ((args*)input)->mat_backup[z+k-1][j]);
		}
	}
	return(0);
 }

void relax(args& Comp_args)
{	
	int size = Comp_args.size;
	int numThreads = Comp_args.numThreads;
	pthread_t* threads = initThreads(numThreads);
	for (int z = 1; z < size-1; ++z)
	{	
		Comp_args.zeile = z;
		pthread_create(&threads[z%numThreads], NULL, &compute, (void *) &Comp_args);
		pthread_join(threads[z%numThreads], NULL);
	}
}

int main(int argc, char* argv[])
{
	//terminal parameter
	int grid_size_N = 1024 + 2; //strtol(argv[1], NULL, 10)+2;	
	int diameter_m = 128; //strtol(argv[2], NULL, 10);
	double initTemp = 100; //strtof(argv[1], NULL);
	int numIter = strtol(argv[1], NULL, 10);
	double initTempCirlce = 50; //strtof(argv[5], NULL);
	int numThreads = strtol(argv[2], NULL, 10);

	double** mat = allocMat(grid_size_N);
	double** mat_backup = allocMat(grid_size_N);
	initMat(mat, grid_size_N, initTemp);
	initCircle(mat,grid_size_N, diameter_m, initTempCirlce);
	copyMatrix(mat, mat_backup, grid_size_N);

	double G[3] = {6./25., -12./25., 6./25.}; // Gewichtungen
	args Comp_args;
	Comp_args.zeile = 0;
	Comp_args.spalte = 0;
	Comp_args.mat = mat;
	Comp_args.mat_backup = mat_backup;
	Comp_args.size = grid_size_N;
	Comp_args.G = G;
	Comp_args.numThreads = numThreads;
	time_t now;
	time_t not_now;
	time_t cylces = 0;
	time_t cylces_ave;
	MPI_Init(NULL, NULL);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	int batch_size = num_iter/world_size;
	if (batch_size % 2)
	{
		printf("only even batch size allowed!\n");
		exit(0);
	}
	 MPI_Status status;
 
 
	for (int itr=0; itr<numIter; ++itr)
	{
		rdtsc(now);
		copyMatrix(mat, mat_backup, grid_size_N);
		relax(Comp_args);
		rdtsc(not_now);
		time_t cycles_new = not_now - now;
		cylces += cycles_new;
		cylces_ave = cylces/(itr + 1);
		printf ("Cylcles fuer diese Iteration  = %ld\n", cycles_new);
		printf ("Cycles Durchschnitt  = %ld\n", cylces_ave);
		printf ("aktuell iter num  = %d\n", itr);
		printf ("iter gesamt  = %d\n", numIter);
		printf ("threads  = %d\n\n", numThreads);
	}
	return(0);
}

