#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <mpi.h>
#include "displaymatrix-mpi.h"
#include "partdiff.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  L               /* # of Lines to generate and calculate           */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

// Zeitmesse
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */

//Global MPI-Data
int rank, numThreads;

static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, double *thread_options)
{
	arguments->N = (thread_options[0] * 8) + 9 - 1;
	arguments->h = 1.0 / arguments->N;
	arguments->L = (rank < arguments->N & numThreads) (arguments->N + numThreads - 1) / numThreads + 1 : (arguments->N + numThreads - 1) / numThreads;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}


// Generiert eine N+1 * L+2 Matrix (Normale Breite, aber aufgeteilte Höhe + Buffer)
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	//TODO: Sichergehen, dass ich hier alles richtig gemacht habe
	// +1 bzw. +2 da es dann unten leserlicher ist.
	uint64_t const N = arguments->N + 1;
	uint64_t const L = arguments->L + 2;

	arguments->M = allocateMemory(2 * N * L * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < 2; i++)
	{
		arguments->Matrix[i] = allocateMemory(N * sizeof(double*));
		
		for (j = 0; j < L; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * N * L) + (j * L);
		}
	}
}

static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < 2; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

static
void
initMatrices (struct calculation_arguments* arguments, double inf_func)
{
	uint64_t g, i, j, limit;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	uint64_t const L = arguments->L;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < 2; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	//TODO: Sichergehen, dass ich hier alles richtig gemacht habe
	if (inf_func == FUNC_F0)
	{
		limit = (rank == numThreads - 1) ? L - 1 : L;
		for (g = 0; g < 2; g++)
		{
			i = (rank == 0) ? 0 : 1;
			for (i; i <= limit; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				if(rank == numThreads - 1) Matrix[g][L][i] = h * i;
				if(rank == 0) Matrix[g][0][i] = 1.0 - (h * i);
			}

			if(rank == numThreads - 1) Matrix[g][L][0] = 0.0;
			if(rank == 0) Matrix[g][0][N] = 0.0;
		}
	}
}

static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, double *thread_options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (L + 1) * sizeof(double) * 2 / 1024.0 / 1024.0);
	printf("Berechnungsmethode: Jacobi\n");

	printf("Interlines:         %" PRIu64 "\n", (int) thread_options[0]);
	printf("Stoerfunktion:      ");

	if (thread_options[1] == FUNC_F0)
		printf("f(x,y) = 0\n");
	else
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)\n");

	printf("Terminierung:       ");

	if (thread_options[2] == TERM_PREC)
		printf("Hinreichende Genaugkeit\n");
	else
		printf("Anzahl der Iterationen\n");

	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

void calculate (struct calculation_arguments *arguments, struct calculation_results *results, double *thread_options)
{
	int target = rank + 1;
	int source = rank - 1;
	int i, j, *done = 0;
	int m1 = 0, m2 = 1;

	int const N = arguments->N;
	int const L = arguments->L;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = thread_options[3];

	double **Matrix_Out = arguments->Matrix[m1];
	double **Matrix_In = arguments->Matrix[m2];
	double *maxresiduum, star = 0.0, residuum = 0.0;
	MPI_Request reqUpper, reqLower;
	
	if(rank == 0 && thread_options[2] == TERM_PREC)
		double *maxres = malloc(sizeof(double));
	else if(thread_options[2] == TERM_ITER)
		int iteration;
	
	while(!*done)
	{
		/* over all rows */
		for (i = 1; i < L - 1; i++)
		{
			if (thread_options[1] == FUNC_FPISIN)
			{
				fpisin_i = *param->fpisin * sin(*param->pih * (double)i);
			}
			//* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
				if (thread_options[1] == FUNC_FPISIN)
				{
					star += fpisin_i * sin(*param->pih * (double)j);
				}
				if (thread_options[2] == TERM_PREC || thread_options[3] == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					*maxresiduum = (residuum < *maxresiduum) ? *maxresiduum : residuum;
				}
				Matrix_Out[i][j] = star;
			}
		}
		if(thread_options[2] == TERM_PREC)
		{
			MPI_Reduce(&maxresiduum, &maxres, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			if(rank == 0)
			{
				if(maxres < thread_options[4])
					done = 1;
			}
			MPI_Bcast(&done, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
		else
		{
			iteration++;
			if(iteration == thread_options[3])
				done = 1;
		}
		// First and last Thread don't send/recv first or last row
		if(rank != 0)
			MPI_Isend(&Matrix_Out[0], L, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &reqUpper);
		if(rank != numThreads - 1)
			MPI_Isend(&Matrix_Out[L-1], L, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &reqLower);
		if(rank != 0)
		{
			MPI_Recv(&Matrix_Out[0], L, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);
			MPI_Wait(&reqUpper);
		}
		if(rank != numThreads - 1)
		{
			MPI_Recv(&Matrix_Out[L-1], L, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);
			MPI_Wait(&reqLower);
		}
		results->stat_iteration++;
		results->stat_precision = maxresiduum;
				
		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;
		Matrix_Out = arguments->Matrix[m1];
		Matrix_In = arguments->Matrix[m2];
	}
	results->m = m2;
}

int main (int argc, char *argv)
{
	MPI_Init(&argc, &argv);
	
	double thread_options[numThreads];
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);

	//make most options double, for better MPI-Sending :)
	if(rank == 0)
	{
		askParams(&options, argc, argv);
		thread_options[0] = (double) options.interlines;
		thread_options[1] = (double) options.inf_func;
		thread_options[2] = (double) options.termination;
		thread_options[3] = (double) options.term_iteration;
		thread_options[4] = (double) options.term_precision;
	}
	MPI_Bcast(&thread_options, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	initVariables(&arguments, &results, &thread_options);
	
	allocateMatrices(&arguments);
	initMatrices(&arguments, thread_options[1]);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &thread_options);
	gettimeofday(&comp_time, NULL);

	MPI_Barrier(MPI_COMM_WORLD);
	
	displayStatistics(&arguments, &results, &thread_options);
	DisplayMatrix (&arguments, &results, &thread_options, rank, (int) arguments->L, 1, (int) arguments->L-1);

	freeMatrices(&arguments);
	
	MPI_Finalize();
	return 0;
}
