#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

//these are some very handsome travellers
int rank, numThreads, length, rest, first;

int*
init (void)
{
	int* buf = malloc(sizeof(int) * length);

	//Niemand darf jemals das gesamte Array im Speicher haben!
	if(rank == 0)
	{
		int i, j;
		srand(time(NULL));
		
		int sendbuf[length]; 
		for (i = 0; i < numThreads; i++)
		{
			for (j = 0; j < length; j++)
			{
				//wenn Reste vorhanden, soll dieser auf die ersten x Threads verteilt werden, nachfolgende Threads haben im letzten Arrayeintrag 'NULL'
				if(i >= numThreads - (numThreads - rest)  && j == length - 1)
					sendbuf[j] = -1;
				else
				{
					// Do not modify "% 25" - we didn't :)
					sendbuf[j] = rand() % 25;
				}
			}

			//Send data, or save, if you are Zero
			if(i != 0)
				MPI_Send(sendbuf, length, MPI_INT, i, 0, MPI_COMM_WORLD);
			else
			{
				for(j = 0; j < length; j++)
					buf[j] = sendbuf[j];
			}
		}

		//Send first int to last Thread
		MPI_Send(&buf[0], 1, MPI_INT, numThreads - 1, 0, MPI_COMM_WORLD);
	}
	else
	{
		//Receive your part of the pie :)
		MPI_Recv(buf, length, MPI_INT, 0, 0, MPI_COMM_WORLD , MPI_STATUS_IGNORE);
	}

	//Last Thread receives the first int.
	if(rank == numThreads - 1)
		MPI_Recv(&first, 1, MPI_INT, 0, 0, MPI_COMM_WORLD , MPI_STATUS_IGNORE);
	
	return buf;
}

void
circle (int* buf)
{
	MPI_Status stat;
	int sendbuf[length];
	int done = 0;
	int target = ((rank + 1) % numThreads);
	int source = (((rank - 1) % numThreads) + numThreads) % numThreads;
	MPI_Request req;
	int i;

	while(done == 0)
	{
		for(i = 0; i < length; i++)
			sendbuf[i] = buf[i];
		//send data asynchronisly, for immediate receive
		MPI_Isend(sendbuf, length, MPI_INT, target, 0, MPI_COMM_WORLD, &req);
		MPI_Recv(buf, length, MPI_INT, source, 0, MPI_COMM_WORLD, &stat);
		
		//wait for sent data and write recvbuf to buf
		MPI_Wait(&req, &stat);

		//if we are done, broadcast so.
		if(rank == numThreads - 1 && buf[0] == first)
			done = 1;
		MPI_Bcast(&done, 1, MPI_INT, numThreads - 1, MPI_COMM_WORLD);
	}
}

int
main (int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int N, i, j;
	int *buf;

	if (argc < 2)
	{
		printf("Arguments error!\nPlease specify a buffer size.\n");
		return EXIT_FAILURE;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);

	// Array Length
	N = atoi(argv[1]);
	length = (N + numThreads - 1) / numThreads;
	rest = N % numThreads;
	// init part of array
	buf = init();
	
	//some printy boys
	if(rank == 0)
		printf("\nBEFORE\n");
	for (i = 0; i < numThreads; i++)
	{	
		if (i == rank)
		{
			for (j = 0; j < length; j++)
			{
				if(buf[j] > -1)
					printf("rank %d: %d\n", rank, buf[j]);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	circle(buf);
	
	//more printy boys :o
	if(rank == 0)
		printf("\nAFTER\n");
	for (i = 0; i < numThreads; i++)
	{	
		if (i == rank)
		{
			for (j = 0; j < length; j++)
			{
				if(buf[j] > -1)
					printf("rank %d: %d\n", rank, buf[j]);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	free(buf);
	MPI_Finalize();
	return EXIT_SUCCESS;
}
