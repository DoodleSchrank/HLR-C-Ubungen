/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
//#include </usr/local/Cellar/mpich/3.3.1/include/mpi.h>

#include "partdiff.h"

struct calculation_arguments {
	uint64_t  N;                        /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;             /* number of matrices                             */
	double    h;                        /* length of a space between two lines            */
	double    ***Matrix;                /* index matrix used for addressing M             */
	double    *M;                       /* two matrices with real values                  */
	uint64_t chunk_size;                /* Größe des zu bearbeitenden Teilabschnitts von jedem Prozess */
	uint64_t chunk_start_index;         /* Startindex des chunks relativ zur Gesamtmatrix */
	uint64_t chunk_end_index;           /* Endindex des chunks relativ zur Gesamtmatrix*/
	uint64_t actual_start_index;        /* Eigentlicher Startindex unter Berücksichtung der benötigten Nachbarzellen */
};

struct calculation_results {
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static void initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options, int rank, int nprocs) {
    uint64_t remainder;

    arguments->N = (options->interlines * 8) + 9 - 1;
    arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
    arguments->h = 1.0 / arguments->N;

    arguments->chunk_size = (arguments->N - 1) / nprocs;                                              /* Berechnet die Chunkgröße */
    remainder = (arguments->N - 1) % nprocs;                                           /* Berechnet den Rest, falls N nicht gleichmäßig auf die Anzahl der Prozesse verteilt werden kann */

    if (remainder != 0 && rank < (int) remainder) {                                             /* Erhöht die Chunkgröße der Prozesse mit Rang < Rest um 1 */
        arguments->chunk_size++;
        arguments->chunk_start_index = arguments->chunk_size * rank + 1;                            /* Globaler Startindex des Chunks mit Rest (ohne Nachbarzellen) */
    } else {
        arguments->chunk_start_index = arguments->chunk_size * rank + remainder + 1;                /* Globaler Startindex des Chunks ohne Rest (ohne Nachbarzellen) */
    }

    arguments->chunk_end_index = arguments->chunk_start_index + arguments->chunk_size - 1;      /* Berechnung des Endindexes */


    arguments->chunk_size += 2;                                                             /* Restlichen Prozesse benötigen zwei Zeilen, jeweils eine von Prozess rank - 1 und eine von Prozess rank + 1;*/


    arguments->actual_start_index = arguments->chunk_start_index - 1;                           /* Berechnung des eigentlichen Startindexes (relevant, falls Störfunktion 1 gewählt wird) */

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static void freeMatrices (struct calculation_arguments* arguments) {
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)   {
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static void* allocateMemory (size_t size) {
	void *p;

	if ((p = malloc(size)) == NULL) {
		printf("Speicherprobleme! (%zu Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static void allocateMatrices (struct calculation_arguments* arguments) {
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * arguments->chunk_size * (N + 1) * sizeof(double));          /* (N + 1) wurde mit Chunkgröße ersetzt */
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++) {
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j < arguments->chunk_size; j++) {                                                                        /* j <= N wurde mit j < Chunkgröße ersetzt */
			arguments->Matrix[i][j] = arguments->M + (i * arguments->chunk_size * (N + 1)) + (j * (N + 1));                  /* (N + 1) wurde mit Chunkgröße ersetzt */
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static void initMatrices (struct calculation_arguments* arguments, struct options const* options, int rank, int nprocs) {
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++) {
        for (i = 0; i < arguments->chunk_size; i++) {                                                                              /* i <= N wurde mit i < Chunkgröße ersetzt */
			for (j = 0; j <= N; j++) {
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0) {
		for (g = 0; g < arguments->num_matrices; g++) {
			for (i = 0; i < arguments->chunk_size; i++) {
				Matrix[g][i][0] = 1.0 - (h * (uint64_t)(i + arguments->actual_start_index));                                  /* Randwert hängt von der globalen Position der Zeile ab, es wird hier jedoch mit actual_start_index gerechnet, da
				                                                                                                               * entsprechend auch die Nachbarzellen des chunks mit Randwerten initialisiert werden müssen */
				Matrix[g][i][N] = h * (uint64_t)(i + arguments->actual_start_index);
			}

			if (rank == 0) {                                                                                                  /* Initialisierung der oberen Randzeile nur bei Prozess 0 notwendig */
                for (i = 0; i <= N; i++) {
                    Matrix[g][0][i] = 1.0 - (h * i);
                }
                Matrix[g][0][N] = 0.0;
            }

			if (rank == (nprocs - 1)) {                                                                                       /* Initialisierung der unteren Randzeile nur bei Prozess nprocs - 1 notwendig */
                for (i = 0; i <= N; i++) {
                    Matrix[g][arguments->chunk_size - 1][i] = h * i;
                }
                Matrix[g][arguments->chunk_size - 1][0] = 0.0;
            }
		}
	}
}

/* ************************************************************************ */
/* calculateSeq: solves the equation                                        */
/* ************************************************************************ */
static void calculateSeq (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options) {
    int i, j;                                   /* local variables for loops */
    int m1, m2;                                 /* used as indices for old and new matrices */
    double star;                                /* four times center value minus 4 neigh.b values */
    double residuum;                            /* residuum of current iteration */
    double maxresiduum;                         /* maximum residuum value of a slave in iteration */

    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

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

    while (term_iteration > 0) {
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In  = arguments->Matrix[m2];

        maxresiduum = 0;

        /* over all rows */
        for (i = 1; i < N; i++) {
            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN) {
                fpisin_i = fpisin * sin(pih * (double)i);
            }

            /* over all columns */
            for (j = 1; j < N; j++) {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

                if (options->inf_func == FUNC_FPISIN) {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1) {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }

                Matrix_Out[i][j] = star;
            }
        }

        results->stat_iteration++;
        results->stat_precision = maxresiduum;

        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;

        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC) {
            if (maxresiduum < options->term_precision) {
                term_iteration = 0;
            }
        } else if (options->termination == TERM_ITER) {
            term_iteration--;
        }
    }

    results->m = m2;
}


/* ************************************************************************ */
/* calculateJacobi: solves the equation                                     */
/* ************************************************************************ */
static void calculateJacobi (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, int rank, int nprocs) {
	int i, j;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */
	double res_buf;                                                         /* Residuum buffer zur Zwischenspeicherung für MPI_Allreduce */
	MPI_Request *requests = malloc(sizeof(MPI_Request) * 4);                /* Request array für MPI_Waitall */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

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

	while (term_iteration > 0) {
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

		/* over all rows */
		for (i = 1; i < (int) arguments->chunk_size - 1; i++) {                                                       /* i < N wurde mit i < chunk_size ersetzt */
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN) {
                    fpisin_i = fpisin * sin(pih * (double)(i + arguments->chunk_start_index - 1));              /*fpisin_i hängt von der globalen Zeilenposition ab, daher addition des globalen Startindexes zu i*/
			}

			/* over all columns */
			for (j = 1; j < N; j++)	{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN) {
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1) {
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}
		if (nprocs > 1) {                                                                                                                   /* Versenden von Nachrichten bei Ausführung mit einem Prozess nicht notwendig */
		    if (rank != 0) {                                                                                                                /* Rang 0 sendet nur seine letzte Zeile und erhält die erste Zeile von Prozess 1 */
                MPI_Isend(Matrix_Out[1], N, MPI_DOUBLE, rank - 1, 10, MPI_COMM_WORLD, &requests[0]);
                MPI_Irecv(Matrix_Out[0], N, MPI_DOUBLE, rank - 1, 20, MPI_COMM_WORLD, &requests[1]);
            }
		    if (rank != (nprocs - 1)) {                                                                                                     /* Rang nprocs - 1 sendet nur seine erste Zeile und erhält die letzte Zeile von Prozess nprocs - 2 */
                MPI_Isend(Matrix_Out[arguments->chunk_size - 2], N, MPI_DOUBLE, rank + 1, 20, MPI_COMM_WORLD, &requests[2]);
                MPI_Irecv(Matrix_Out[arguments->chunk_size - 1], N, MPI_DOUBLE, rank + 1, 10, MPI_COMM_WORLD, &requests[3]);
            }

		    if (rank == 0) {
		        MPI_Waitall(2, &requests[2], MPI_STATUSES_IGNORE);                                                   /* Rang 0 muss nur auf die ersten beiden requests warten */
            } else if (rank == (nprocs - 1)) {
                MPI_Waitall(2, &requests[0], MPI_STATUSES_IGNORE);                                                   /* Rang nprocs - 1 muss nur auf die letzten beiden requests warten */
            } else {
                MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);                                                                       /* Alle anderen Prozesse müssen auf alle vier requests warten */
            }
        }

		MPI_Allreduce(&maxresiduum, &res_buf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);                                               /* Ermittlung von maxresiduum von allen Prozessen und Übermittlung des Ergebnisses an alle Prozesse */
		maxresiduum = res_buf;                                                                                                       /* Belegung von maxresiduum mit dem globalen maxresiduum */

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC) {
			if (maxresiduum < options->term_precision) {
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER) {
			term_iteration--;
		}
	}

	free(requests);                                                                                                                         /* Löschen des requests Arrays */
	results->m = m2;
}

/* ************************************************************************ */
/* calculateGaussSeidel: solves the equation                                 */
/* ************************************************************************ */
static void calculateGaussSeidel (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options, int rank, int nprocs) {
    int i, j;                                                       /* local variables for loops */
    int m1 = 0, m2 = 0;                                             /* used as indices for old and new matrices */
    double star;                                                    /* four times center value minus 4 neigh.b values */
    double residuum;                                                /* residuum of current iteration */
    double maxresiduum;                                             /* maximum residuum value of a slave in iteration */

    int const N = arguments->N;
    double const h = arguments->h;

    double pih = 0.0;
    double fpisin = 0.0;

    int term_iteration = options->term_iteration;

    if (options->inf_func == FUNC_FPISIN) {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    // TODO BEGIN (added):
    int finalized = 0 == 1;                                        /* indicates when to finalize calculation */
    MPI_Request request;
    char uselessBuffer;
    if(rank == 0)
        MPI_Irecv(&uselessBuffer, 1, MPI_BYTE, nprocs - 1, 3, MPI_COMM_WORLD, &request);
    // TODO END (added)

    while (!finalized) { // TODO changed line
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In  = arguments->Matrix[m2];

        maxresiduum = 0;
        // TODO BEGIN (added)
        if (rank > 0) {
            MPI_Send(Matrix_In[1], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(Matrix_In[0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&maxresiduum, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //printf("[%i] Startschuss für nächste Iteration!\n", rank);
            finalized = (maxresiduum == -1.0);
        }
        // TODO END (added)

        /* over all rows */
        for (i = 1; i < (int) arguments->chunk_size - 1; i++) {                                                       /* Iteration bis zur vorletzten Speicherzeile */

            // TODO BEGIN (added)
            if (i == (int) arguments->chunk_size - 2)        // Ausführung vor Berechnung der letzten Zeile
                if (rank < nprocs - 1)
                    MPI_Recv(Matrix_In[arguments->chunk_size - 1], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // TODO END (added)


            double fpisin_i = 0.0;

            if (options->inf_func == FUNC_FPISIN) {
                fpisin_i = fpisin * sin(pih * (double)(i + arguments->chunk_start_index - 1));              /*fpisin_i hängt von der globalen Zeilenposition ab, daher addition des globalen Startindexes zu i*/
            }

            /* over all columns */
            for (j = 1; j < N; j++) {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

                if (options->inf_func == FUNC_FPISIN) {
                    star += fpisin_i * sin(pih * (double)j);
                }

                if (options->termination == TERM_PREC || term_iteration == 1) {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }

                Matrix_Out[i][j] = star;
            }
        }

        results->stat_iteration++;
        results->stat_precision = maxresiduum;

        /* exchange m1 and m2 */
        //i = m1;
        //m1 = m2;
        //m2 = i;
        // TODO BEGIN (added)
        if (rank == 0) {
            MPI_Test(&request, &finalized, MPI_STATUS_IGNORE);
        }
        if (rank < nprocs - 1) {
            MPI_Send(Matrix_In[arguments->chunk_size - 2], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            maxresiduum = finalized ? -1.0: maxresiduum;
            if(maxresiduum == -1.0)
                printf("Termination forwarding...!\n");
            MPI_Send(&maxresiduum, 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
        }
        if(rank == nprocs - 1 && options->termination == TERM_PREC && maxresiduum < options->term_precision) {
            MPI_Isend(&uselessBuffer, 1, MPI_BYTE, 0, 3, MPI_COMM_WORLD, &request);
            printf("Termination initiated!\n");
        }
        // TODO END (added)

        /* check for stopping calculation depending on termination method */
        if (options->termination == TERM_PREC) {
            //printf("Not supported!\n");
            //exit(-88);
            //if (maxresiduum < options->term_precision) {
            //    term_iteration = 0;
            //}
        } else if (options->termination == TERM_ITER) {
            term_iteration--;
            finalized = (term_iteration <= 0);
        }
    }
                                                                                                                      /* Löschen des requests Arrays */
    results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static void displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options) {
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL) {
		printf("Gauß-Seidel");
	} else if (options->method == METH_JACOBI) {
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0) {
		printf("f(x,y) = 0");
	} else if (options->inf_func == FUNC_FPISIN) {
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC) {
		printf("Hinreichende Genaugkeit");
	} else if (options->termination == TERM_ITER) {
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static void displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to) {
    int const elements = 8 * options->interlines + 9;

    int x, y;
    double** Matrix = arguments->Matrix[results->m];
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
            /* check whether this line belongs to rank 0 */
            if (line < from || line > to) {
                /* use the tag to receive the lines in the correct order
                 * the line is stored in Matrix[0], because we do not need it anymore */
                MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
            }
        } else {
            if (line >= from && line <= to) {
                /* if the line belongs to this process, send it to rank 0
                 * (line - from + 1) is used to calculate the correct local address */
                MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
            }
        }

        if (rank == 0) {
            for (x = 0; x < 9; x++) {
                int col = x * (options->interlines + 1);

                if (line >= from && line <= to) {
                    /* this line belongs to rank 0 */
                    printf("%7.4f", Matrix[line][col]);
                } else {
                    /* this line belongs to another rank and was received above */
                    printf("%7.4f", Matrix[0][col]);
                }
            }

            printf("\n");
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    fflush(stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int main (int argc, char** argv) {
    struct options options;
    struct calculation_arguments arguments;
    struct calculation_results results;

    int rank = 0;
    int nprocs = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    askParams(&options, argc, argv, rank);

    if((int) ((options.interlines * 8) + 9 - 1) <= nprocs) {
        nprocs = 1;
    }

    initVariables(&arguments, &results, &options, rank, nprocs);
    allocateMatrices(&arguments);
    initMatrices(&arguments, &options, rank, nprocs);

    if (rank == 0)                                                      /* Prozess 0 misst die Berechnungszeit */
        gettimeofday(&start_time, NULL);

    MPI_Barrier(MPI_COMM_WORLD);                                        /* Prozesse warten bis Prozess 0 die Zeitmessung gestartet hat*/

    if (nprocs == 1) {
        if (rank == 0)
            calculateSeq(&arguments, &results, &options);
    } else if (options.method == METH_JACOBI) {
        calculateJacobi(&arguments, &results, &options, rank, nprocs);
    } else {
        calculateGaussSeidel(&arguments, &results, &options, rank, nprocs);
    }
	MPI_Barrier(MPI_COMM_WORLD);                                        /* Barrier bis alle Prozesse mit der Berechnung fertig sind*/
	if (rank == 0) {                                                    /* Zeit wird von Prozess 0 gestoppt und Statistik wird ausgegeben */
        gettimeofday(&comp_time, NULL);
        displayStatistics(&arguments, &results, &options);
    }

    displayMatrix(&arguments, &results, &options, rank, nprocs, arguments.chunk_start_index, arguments.chunk_end_index);

	freeMatrices(&arguments);

	return 0;
}
