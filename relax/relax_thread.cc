#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"
#include <pthread.h>

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double** initMat(int size, double init_val){
	int i, j; 
  	// Allocation of Matrix
    double **mat = (double **)malloc(size * sizeof(double *)); 
    for (i=0; i<size; i++)
    	mat[i] = (double *)malloc(size * sizeof(double)); 
    // Note that arr[i][j] is same as *(*(arr+i)+j) 
    
    for (i = 0; i <  size; i++)
    {
		for (j = 0; j < size; j++)
		{
			if(j+1==size || i+1 == size || j== 0 || i == 0)		mat[i][j] = 0;
			else	mat[i][j] = init_val;
		}
	}
   return mat; 
}

double** init_temp_mat(int size){
	int i, j; 
    double **mat = (double **)malloc(size * sizeof(double *)); 
    for (i=0; i<size; i++)
    	mat[i] = (double *)malloc(size * sizeof(double)); 
    
    for (i = 0; i <  size; i++)
    {
		for (j = 0; j < size; j++)
		{
			mat[i][j] = 0;
		}
	}
   return mat; 
}

void init_circle(double** mat, int param_grid_size, int diameter_m, double H)
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

void copyMatrix(double** mat, double** temp_mat , int size)
{	
	for(int z=0; z<size; ++z)
	{
		for(int s=0; s<size; ++s)
			temp_mat[z][s] = mat[z][s];
	}
}	



pthread_t* initThreads(int count_t)
{
	pthread_t* threads =  (pthread_t* )malloc(count_t * sizeof(pthread_t));
	return(threads);
}

//struct für die Parameterübergabe an Threads
struct args {
		int zeile;
		int spalte;
		double** mat;
		double** temp_mat;
		int size;
		double* G;
		
	};
/*die paralelisierungsmethode für den Thradaufruf
 * es wird für jeden Abstand das Integral berechnet 
 * und auf das gesamt ergebniss addiert
 */
 
 void *compute(void *input){
	int z = ((struct args*)input)->zeile;
	z= z-5;
	double** temp_mat = ((struct args*)input)->temp_mat;
	int size = ((struct args*)input)->size;
	//((struct args*)input)->G;
	for(int j = 1; j < size -5 ; ++j){
	for(int k = 0; k < 3; ++k){
	((struct args*)input)->mat[z][j] +=((struct args*)input)->G[k]*(temp_mat[z][j+k-1] + temp_mat[z+k-1][j]);
	}
	}
	return(0);
 }

void relax(double** mat, double** temp_mat , int size, int count_t)
{	
	copyMatrix(mat, temp_mat, size);
	
	//struct aus Aurgumenten für den thread-aufruf
	struct args *Comp_args = (struct args *)malloc(sizeof(struct args));
	Comp_args->mat = mat;
	Comp_args->temp_mat = temp_mat;
	Comp_args->size = size;
	
	Comp_args-> G[0] = 6/25;
	Comp_args-> G[1] = -12/25;
	Comp_args-> G[2] = 6/25;
	pthread_t* threads = initThreads(count_t);
	for (int z = 1; z < size-1; ++z)
	{	
		Comp_args->zeile = z;
		pthread_create(&threads[z%count_t], NULL, &compute, (void *) Comp_args);
		pthread_join(threads[z%count_t], NULL);
	}
}

int main(int argc, char* argv[])
{
	printf("argc %d", argc);
	//terminal parameter
	int grid_size_N = 1024; //strtol(argv[1], NULL, 10)+2;	
	int diameter_m = 128; //strtol(argv[2], NULL, 10);
	double H = 100; //strtof(argv[1], NULL);
	int iter = strtol(argv[1], NULL, 10);
	double init_val_H = 0; //strtof(argv[5], NULL);
	double count_t = strtof(argv[2], NULL);
	
	double** mat = initMat(grid_size_N, init_val_H );
	double** temp_mat = init_temp_mat(grid_size_N);
	init_circle(mat,grid_size_N, diameter_m, H);



	time_t now;
	time_t not_now;
	time_t cylces = 0;
	time_t cylces_ave;
	
	for (int itr=0; itr<iter; itr++)
	{
		rdtsc(now);
		relax(mat, temp_mat, grid_size_N, count_t);
		rdtsc(not_now);
		time_t cycles_new = not_now - now;
		cylces += cycles_new;
		cylces_ave = cylces/(itr + 1);
		printf ("Cylcles fuer diese Iteration  = %ld\n", cycles_new);
		printf ("Cycles Durchschnitt  = %ld\n", cylces_ave);
	}


	return(0);
}

