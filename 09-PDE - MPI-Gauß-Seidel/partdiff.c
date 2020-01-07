#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <float.h>
#include <omp.h>
#include <mpi.h>

#include "partdiff.h"

// foobar

struct calculation_arguments {
	uint64_t N; /* number of spaces between lines (lines=N+1)	 */
	uint64_t num_matrices; /* number of matrices							 */
	double h; /* length of a space between two lines			*/
	double ** * Matrix; /* index matrix used for addressing M			 */
	double * M; /* two matrices with real values				  */
};

struct calculation_results {
	uint64_t m;
	uint64_t stat_iteration; /* number of current iteration					*/
	double stat_precision; /* actual precision of all slaves in iteration	*/
};

/* ************************************************************************ */
/* Global variables														 */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started					  */
struct timeval comp_time; /* time when calculation completed				*/

//Global MPI-Data
int rank, numThreads;

//Größe und Anfang & Ende der Teilmatrizen der einzelnen Prozesse
uint64_t matrix_size, matrix_from, matrix_to;

/* ************************************************************************ */
/* initVariables: Initializes some global variables						 */
/* ************************************************************************ */
static
void
initVariables(struct calculation_arguments * arguments, struct calculation_results * results, struct options
	const * options) {
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;

	uint64_t N = arguments->N;

	// safety if interlines > numThreads / 4 
	// division by 4 to guarantee at least some benefit from paralellization
	numThreads = (numThreads > (N - 1) / 4) ? numThreads : floor((N - 1) / 4);
	matrix_size = ceil((float)(N - 1) / numThreads);
	matrix_from = ((uint64_t)(matrix_size * rank + 1) < N) ? matrix_size * rank + 1 : N;
	matrix_to = ((uint64_t)(matrix_size * (rank + 1)) < (N - 1)) ? matrix_size * (rank + 1) : N - 1;
	if (matrix_from > matrix_to) {
		matrix_from = N;
		matrix_to = N - 1;
	}
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices								  */
/* ************************************************************************ */
static
void
freeMatrices(struct calculation_arguments * arguments) {
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++) {
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()														*/
/* allocates memory and quits if there was a memory allocation problem	  */
/* ************************************************************************ */
static
void *
	allocateMemory(size_t size) {
		void * p;

		if ((p = malloc(size)) == NULL) {
			printf("Speicherprobleme! (%"
				PRIu64 " Bytes angefordert)\n", size);
			exit(1);
		}

		return p;
	}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices						  */
/* ************************************************************************ */
static
void
allocateMatrices(struct calculation_arguments * arguments) {
	uint64_t i, j;

	uint64_t
	const N = arguments->N;
	uint64_t
	const size = matrix_size + 2;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * size * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double ** ));

	for (i = 0; i < arguments->num_matrices; i++) {
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double * ));
		for (j = 0; j <= N; j++)
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * size) + (j * (N + 1));
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables	   */
/* ************************************************************************ */
static
void
initMatrices(struct calculation_arguments * arguments, struct options const *options) {
	// it fit many lööp, brøther
	uint64_t g, i, j;

	uint64_t const N = arguments->N;
	uint64_t const size = matrix_size + 2;

	double const h = arguments->h;
	double ** * Matrix = arguments->Matrix;

	// initialize matrix/matrices with zeros
	for (g = 0; g < arguments->num_matrices; g++) {
		for (i = 0; i < size; i++) {
			for (j = 0; j < N; j++) {
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	// initialize borders, depending on function (function 2: nothing to do)
	if (options->inf_func == FUNC_F0) {
		for (g = 0; g < arguments->num_matrices; g++) {
			for (i = 0; i < size; i++) {
				Matrix[g][i][0] = 1.0 - (h * (i + matrix_from - 1));
				Matrix[g][i][N] = h * (i + matrix_from - 1);
			}

			if (matrix_from == 1) {
				for (j = 0; j <= N; j++)
					Matrix[g][0][j] = 1.0 - (h * j);
			}

			if (matrix_to >= (N - 1)) {
				for (j = 0; j < N; j++)
					Matrix[g][matrix_size + 1][j] = h * j;
			}
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation										   */
/* ************************************************************************ */
static
void
calculate(struct calculation_arguments
	const * arguments, struct calculation_results * results, struct options
	const * options) {
	int target = rank + 1;
	int source = rank - 1;
	uint64_t i, j = 0;
	int m1, m2;

	uint64_t
	const N = arguments->N;
	double
	const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	double maxresiduum, star, residuum;
	double maxres = DBL_MAX;
	MPI_Request reqSendFirst, reqSendLast, reqRecvFirst, reqRecvLast, reqRes;

	double ** Matrix_Out;
	double ** Matrix_In;

	double maxresidaa[numThreads];

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI) {
		m1 = 0;
		m2 = 1;
	} else {
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN) {
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	//rank 0 prepares for receiving maxresidaä of other processes
	if (rank == 0 && options->termination == TERM_PREC) {
		if (numThreads > 1) {
			for (i = 1; i < numThreads; i++)
				MPI_Irecv( & maxresidaa[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, & reqRes);
		}
	}
	//the others prepare for receiving maxres (the one that cancels this whole operation)
	else if (options->termination == TERM_PREC)
		MPI_Irecv( & maxres, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, & reqRes);

	// Recieve first row if Gauß-Seidel, needs to be done before lööp
	if (rank != 0 && options->method == METH_GAUSS_SEIDEL)
		MPI_Irecv(Matrix_In[0], N + 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, & reqRecvFirst);

	while (term_iteration > 0) {
		Matrix_Out = arguments->Matrix[m1];
		Matrix_In = arguments->Matrix[m2];

		maxresiduum = 0;

		omp_set_dynamic(0);
		#pragma omp parallel for private(j, star, residuum) reduction(max: maxresiduum) num_threads(options->number)
		for (i = 1; i < N; i++) {
			double fpisin_i = 0.0;
			if (options->inf_func == FUNC_FPISIN)
				fpisin_i = fpisin * sin(pih * (double) i);

			// Wait for first row to be recieved
			if (rank > 0 && i == 0) {
				MPI_Wait( & reqSendFirst, MPI_STATUS_IGNORE);
				MPI_Wait( & reqRecvFirst, MPI_STATUS_IGNORE);
			}
			// First row to be sent
			if (rank > 0 && i == 1) {
				MPI_Isend(Matrix_Out[1], N + 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, & reqSendFirst);
				MPI_Irecv(Matrix_Out[0], N + 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, & reqRecvFirst);
			}

			// Wait for last row to be sent
			if (results->stat_iteration > 0 && rank < numThreads - 1 && i == matrix_size) {
				MPI_Wait( & reqSendLast, MPI_STATUS_IGNORE);
				MPI_Wait( & reqRecvLast, MPI_STATUS_IGNORE);
			}

			for (j = 1; j < N; j++) {
				star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] + Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);

				if (options->inf_func == FUNC_FPISIN)
					star += fpisin_i * sin(pih * (double) j);

				if (options->termination == TERM_PREC || term_iteration == 1) {
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}
		// Send last row
		// ignore last rank because it has no followers /BIG SAD/
		if (rank < numThreads - 1) {
			MPI_Isend(Matrix_Out[matrix_size], N + 1, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, & reqSendLast);
			MPI_Irecv(Matrix_Out[matrix_size + 1], N + 1, MPI_DOUBLE, target, 0, MPI_COMM_WORLD, & reqRecvLast);
		}

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC) {
			//Reduce Maxresiduum and, if criteria is met, send maxres to other processes to kill them (after some time)
			if (rank == 0) {
				maxresidaa[0] = maxresiduum;
				maxres = 0.0;
				for (i = 0; i < numThreads; i++)
					maxres = (maxresidaa[i] > maxres) ? maxresidaa[i] : maxres;
				if (maxres < options->term_precision) {
					for (i = 1; i < numThreads; i++)
						MPI_Isend( & maxres, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, & reqRes);
				}
			} else
				MPI_Isend( & maxresiduum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, & reqRes);
			if (maxres < options->term_precision)
				term_iteration = 0;
		} else if (options->termination == TERM_ITER) {
			term_iteration--;
		}
		results->stat_iteration++;
		results->stat_precision = maxresiduum;
	}

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation	   */
/* ************************************************************************ */
static
void
displayStatistics(struct calculation_arguments
	const * arguments, struct calculation_results
	const * results, struct options
	const * options) {
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:	%f s \n", time);
	printf("Speicherbedarf:	 %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL) {
		printf("Gauß-Seidel");
	} else if (options->method == METH_JACOBI) {
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:		 %"
		PRIu64 "\n", options->interlines);
	printf("Stoerfunktion:	  ");

	if (options->inf_func == FUNC_F0) {
		printf("f(x,y) = 0");
	} else if (options->inf_func == FUNC_FPISIN) {
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:	   ");

	if (options->termination == TERM_PREC) {
		printf("Hinreichende Genaugkeit");
	} else if (options->termination == TERM_ITER) {
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %"
		PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:							   **/
/**																		**/
/** Die Funktion displayMatrix gibt eine Matrix							**/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**																		**/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix	**/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.									  **/
/****************************************************************************/
static
void
DisplayMatrix(struct calculation_arguments * arguments, struct calculation_results * results, struct options * options, int rank, int size, int from, int to) {
	int
	const elements = 8 * options->interlines + 9;

	int x, y;
	double ** Matrix = arguments->Matrix[results->m];
	MPI_Status status;

	/* first line belongs to rank 0 */
	if (rank == 0)
		from--;

	/* last line belongs to rank size - 1 */
	if (rank + 1 == size)
		to++;

	if (rank == 0)
		printf("Matrix:\n");

	for (y = 0; y < 9; y++) {
		int line = y * (options->interlines + 1);
		if (rank == 0) {
			if (line < from || line > to)
				MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, & status);
		} else {
			if (line >= from && line <= to)
				MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
		}

		if (rank == 0) {
			for (x = 0; x < 9; x++) {
				int col = x * (options->interlines + 1);
				if (line >= from && line <= to)
					printf("%7.4f", Matrix[line][col]);
				else
					printf("%7.4f", Matrix[0][col]);
			}
			printf("\n");
		}
	}
	fflush(stdout);
}

/* ************************************************************************ */
/*  main																	*/
/* ************************************************************************ */
int
main(int argc, char ** argv) {
	MPI_Init( & argc, & argv);
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	MPI_Comm_rank(MPI_COMM_WORLD, & rank);
	MPI_Comm_size(MPI_COMM_WORLD, & numThreads);
	askParams( & options, argc, argv, rank);

	initVariables( & arguments, & results, & options);

	// if numthreads > interlines, numthreads gets reduced
	// thus some threads dont need to work as much
	//	- lazy bastards
	if (rank <= numThreads) {
		allocateMatrices( & arguments);
		initMatrices( & arguments, & options);
		gettimeofday( & start_time, NULL);
		calculate( & arguments, & results, & options);
		gettimeofday( & comp_time, NULL);
		MPI_Barrier(MPI_COMM_WORLD);

		// Nur Rang 0 gibt die Statistiken aus
		if (rank == 0)
			displayStatistics( & arguments, & results, & options);

		// Für die Matrix muss allerdings jeder was abliefern
		DisplayMatrix( & arguments, & results, & options, rank, numThreads, matrix_from, matrix_to);
		freeMatrices( & arguments);
	}

	MPI_Finalize();
	return 0;
}
