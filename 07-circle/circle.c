#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

int*
init (int length, int rest, int numThreads, int rank, int first)
{
	int* buf = malloc(sizeof(int) * length);

	//Niemand darf jemals das gesamte Array im Speicher haben!
	if(rank == 0)
	{
		int i, j;
		srand(time(NULL));
		
		int* sendbuf = malloc(sizeof(int) * length);
		for (i = 0; i < numThreads; i++)
		{
			for (j = 0; j < length; j++)
			{
				// Do not modify "% 25" - we didn't :)
				sendbuf[i] = rand() % 25;
				//wenn Reste vorhanden, soll dieser auf die ersten x Threads verteilt werden, nachfolgende Threads haben im letzten Arrayeintrag 'NULL'
				if(i >= numThreads - rest && j == length - 1)
					sendbuf[i] == NULL;
			}
			//Send data, or save, if you are Zero
			if(i == 0)
				buf = sendbuf;
			else
				MPI_Send(&buf, length, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		//Send first int to last Thread
		MPI_Send(&buf[0], 1, MPI_INT, numThreads - 1, 0, MPI_COMM_WORLD);
	}
	else
	{
		//Receive your part of the pie :)
		MPI_Recv(&buf, length, MPI_INT, 0, 0, MPI_COMM_WORLD , MPI_STATUS_IGNORE);
	}
	//Last Thread receives the first int.
	if(rank == numThreads - 1)
		MPI_Recv(&first, 1, MPI_INT, 0, 0, MPI_COMM_WORLD , MPI_STATUS_IGNORE);

	return buf;
}

int*
circle (int* buf, int numThreads, int rank, int length, int first)
{
	int* sendbuf = malloc(sizeof(int) * length);
	int done = 0;
	int target = (rank + 1) % numThreads;
	int source = (rank - 1) % numThreads;
	MPI_Request req;
	
	while(done == 0)
	{
		savebuf = buf;
		
		//send data asynchronisly, for immediate receive
		MPI_Isend(&sendbuf, length, MPI_INT, target, 0, MPI_COMM_WORLD, &req);
		MPI_Recv(&buf, length, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//wait for sent data
		MPI_Wait(&req, MPI_STATUS_IGNORE);
		
		//if we are done, broadcast so.
		if(rank == numThreads - 1 && buf[0] == first)
			done = 1;
		MPI_Bcast(&done, 1, MPI_INT, numThreads - 1, MPI_COMM_WORLD);
	}
	return buf;
}

int
main (int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int N, rank, numThreads, length, rest, first, i, j;
	int *buf;

	if (argc < 2)
	{
		printf("Arguments error!\nPlease specify a buffer size.\n");
		return EXIT_FAILURE;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);

	// Array Length
	N = atoi(argv[1]);
	length = N / numThreads + 1;
	rest = N % numThreads;
	buf = init(length, rest, numThreads, rank, first);

	printf("\nBEFORE\n");

	for (i = 0; i < length; i++)
	{
		if(buf[i])
			printf("rank %d: %d\n", rank, buf[i]);
	}

	circle(buf, numThreads, rank, length);

	printf("\nAFTER\n");

	for (i = 0; i < length; i++)
	{
		if(buf[i])
			printf("rank %d: %d\n", rank, buf[i]);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
