#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <mpi.h>
#include "partdiff.h"

struct calculation_arguments
{
	uint64_t  N;			  /* number of spaces between lines (lines=N+1)	 */
	uint64_t  num_matrices;   /* number of matrices							 */
	double	h;			  /* length of a space between two lines			*/
	double	**Matrix;	  /* index matrix used for addressing M			 */
	double	*M;			 /* two matrices with real values				  */
};

struct calculation_results
{
	uint64_t  stat_iteration; /* number of current iteration					*/
	double	stat_precision; /* actual precision of all slaves in iteration	*/
};

// Zeitmesse
struct timeval start_time;	   /* time when program started					  */
struct timeval comp_time;		/* time when calculation completed				*/

//Global MPI-Data
int rank, numThreads;

//Größe und Anfang & Ende der Teilmatrizen der einzelnen Prozesse
uint64_t matrix_size, matrix_from, matrix_to;

static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = 2;
	arguments->h = 1.0 / arguments->N;

	results->stat_iteration = 0;
	results->stat_precision = 0;

	uint64_t N = arguments->N;
	
	//safety if interlines > numThreads / 4 
	//division by 4 to guarantee at least some benefit from paralellization
	numThreads = (numThreads > (N - 1.0) / 4.0) ? numThreads : floor((N-1) / 4);
  matrix_size =  ceil((float)(N-1) / numThreads);
	matrix_from = ((uint64_t) (matrix_size * rank + 1) < N) ? matrix_size * rank + 1 : N;
	matrix_to = ((uint64_t) (matrix_size * (rank + 1)) < (N - 1)) ? matrix_size * (rank + 1) : N - 1;
  if (matrix_from > matrix_to)
  {
	  matrix_from = N;
	  matrix_to = N - 1;
  }
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


// Generiert eine Matrix (Normale Breite, aber aufgeteilte Höhe + Buffer)
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	uint64_t const N = arguments->N;
	uint64_t const size = matrix_size + 2;

	arguments->M = allocateMemory((N + 1) * size * sizeof(double));

	arguments->Matrix = allocateMemory((N + 1) * sizeof(double*));

	for (i = 0; i <= N; i++)
	{
		arguments->Matrix[i] = arguments->M + (i * (N + 1) * size) + (i * (N + 1));
	}
}

static
void
freeMatrices (struct calculation_arguments* arguments)
{
	free(arguments->Matrix);
	free(arguments->M);
}

static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	// it fit many lööp, brøther
	uint64_t i, j;

	uint64_t const N = arguments->N;
	uint64_t const size = matrix_size + 2;

	double const h = arguments->h;
	double** Matrix = arguments->Matrix;

	// initialize matrix/matrices with zeros
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < N; j++)
		{
			Matrix[i][j] = 0.0;
		}
	}
	
	// initialize borders, depending on function (function 2: nothing to do)
	if (options->inf_func == FUNC_F0)
	{
		for (i = 0; i < size; i++)
		{
			Matrix[i][0] = 1.0 - (h * (i + matrix_from - 1));
			Matrix[i][N] = h * (i + matrix_from - 1);
		}
	
		if (matrix_from == 1)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[0][j] = 1.0 - (h * j);
			}
		}
	
		if (matrix_to >= (N - 1))
		{
			for (j = 0; j < N; j++)
			{
				Matrix[matrix_size + 1][j] = h * j;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:	%f s \n", time);
	printf("Speicherbedarf:	 %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: Jacobi\n");

	printf("Interlines:		 %" PRIu64 "\n", options->interlines);
	printf("Stoerfunktion:	  ");

	if (options->inf_func == FUNC_F0)
		printf("f(x,y) = 0\n");
	else
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)\n");

	printf("Terminierung:	   ");

	if (options->termination == TERM_PREC)
		printf("Hinreichende Genaugkeit\n");
	else
		printf("Anzahl der Iterationen\n");

	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

void calculate (struct calculation_arguments *arguments, struct calculation_results *results,  struct options const* options)
{
	int target = 0;//rank + 1;
	int source = 0;//rank - 1;
	uint64_t i, j = 0;

	uint64_t const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	double maxresiduum, star, residuum, *maxres;
	MPI_Status status;
	MPI_Request reqUpper, reqRes;

	double** Matrix = arguments->Matrix;
	double maxresidaa[numThreads];

	double fpisin_i = 0.0;
	
	//rank 0 prepares for receiving maxresidaa of other processes
	if(rank == 0 && numThreads > 1)
	{
		for(i = 1; i < numThreads; i++)
			MPI_Irecv(&maxresidaa[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &reqRes);
	}
	//the others prepare for receiving maxres (the one that cancels this whole operation)
	else
		MPI_Irecv(&maxres, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &reqUpper);
	if (options->termination == TERM_PREC && rank != 0)
	
		MPI_Irecv(&maxres, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &reqRes);
	

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while(term_iteration > 0)
	{
		maxresiduum = 0;
		
		//Only wait after first iteration, else we're stuck
		if(results->stat_iteration > 0 && rank < numThreads - 1)
			MPI_Wait(&reqUpper, MPI_STATUS_IGNORE);
		// ignore rank 0 because it has no previous partner (it's still single and lookin' for a soulmate)
		
		for (i = 1; i < matrix_size + 1; i++)
		{
			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (i + matrix_from - 1));
			}
			
			if (source > 0 && i == matrix_size)
			{
				MPI_Recv(Matrix[0],  N + 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
			}
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix[i-1][j] + Matrix[i][j-1] + Matrix[i][j+1] + Matrix[i+1][j]);
				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * j);
				}
				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}
				Matrix[i][j] = star;
			}
		}
		


		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		// ignore last rank because it has zero followers (on instagram) /BIG SAD/
		if (target < numThreads - 1)
		{
			MPI_Isend(Matrix[matrix_size], N + 1, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, &reqUpper);
		}
		
		if (options->termination == TERM_PREC)
		{	
			//Reduce Maxresiduum and, if criteria is met, send maxres to other processes to kill them (after some time)
			if(rank == 0)
			{
				maxresidaa[0] = maxresiduum;
				for(i = 0; i < numThreads; i++)
				{
					*maxres = (maxresidaa[i] > *maxres) ? maxresidaa[i] : *maxres;
				}
				if(maxresidaa[i] > options->term_precision)
				{
					for(i = 1; i < numThreads; i++)
					{
						MPI_Isend(&maxres, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &reqRes);
					}
				}
			}
			if (*maxres < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}
}

/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
	int const elements = 8 * options->interlines + 9;

	int x, y;
	double** Matrix = arguments->Matrix;
	MPI_Status status;

	/* first line belongs to rank 0 */
	if (rank == 0)
		from--;

	/* last line belongs to rank size - 1 */
	if (rank + 1 == size)
		to++;

	if (rank == 0)
		printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		int line = y * (options->interlines + 1);

		if (rank == 0)
		{
			/* check whether this line belongs to rank 0 */
			if (line < from || line > to)
			{
				/* use the tag to receive the lines in the correct order
				 * the line is stored in Matrix[0], because we do not need it anymore */
				MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
			}
		}
		else
		{
			if (line >= from && line <= to)
			{
				/* if the line belongs to this process, send it to rank 0
				 * (line - from + 1) is used to calculate the correct local address */
				MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
			}
		}

		if (rank == 0)
		{
			for (x = 0; x < 9; x++)
			{
				int col = x * (options->interlines + 1);

				if (line >= from && line <= to)
				{
					/* this line belongs to rank 0 */
					printf("%7.4f", Matrix[line][col]);
				}
				else
				{
					/* this line belongs to another rank and was received above */
					printf("%7.4f", Matrix[0][col]);
				}
			}

			printf("\n");
		}
	}

	fflush(stdout);
}

int main (int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);

	//Rang 0 macht den ganzen Kram
	askParams(&options, argc, argv, rank);

	initVariables(&arguments, &results, &options);
	
	// if numthreads > interlines, numthreads gets reduced
	// thus some threads dont need to work as much
	//	- lazy bastards
	if(rank <= numThreads)
	{
		allocateMatrices(&arguments);
		initMatrices(&arguments, &options);

		gettimeofday(&start_time, NULL);
		calculate(&arguments, &results, &options);
		gettimeofday(&comp_time, NULL);

		MPI_Barrier(MPI_COMM_WORLD);
		
		// Nur Rang 0 gibt die Statistiken aus
		if (rank == 0)
		{
			displayStatistics(&arguments, &results, &options);
		}
		//Für die Matrix muss allerdings jeder was abliefern
		DisplayMatrix (&arguments, &results, &options, rank, numThreads, matrix_from, matrix_to);
	}

	freeMatrices(&arguments);

	MPI_Finalize();
	return 0;
}
