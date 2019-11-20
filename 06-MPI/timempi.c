#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>

void meister(int numThreads)
{
	MPI_Datatype datastruct;
	for(int i = 1; i < numThreads; i++)
	{
		MPI_Recv(&datastruct, 1, datastruct, i, void , void , void);
		printf("%s:%ld\n", datastruct.hostname, datastruct.time);
	}
	
}

void sklave()
{
	MPI_Datatype datastruct;
	char hostname[256];
	struct timeval tv;
	struct timezone tz;
	
	if(gethostname(hostname, sizeof(hostname)) != 0) hostname = "Error";
	if(gettimeofday(&tv, &tz) != 0) tv.tv_usec = 0;
	
	struct data
	{
		char *hostname = hostname;
		long int time = tv.tv_usec;
	} data;
	MPI_Send(&data, 1, datastruct, 0, void, void);
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int myid, numThreads;
	MPI_Datatype datastruct;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
	int blocklength[] = {1, 1};
	int displacements[] = {0, 1};
	MPI_Datatype types[] = {MPI_CHAR, MPI_LONG};

	MPI_Type_struct(2, blocklength, displacements, types, &datastruct);
	MPI_Type_commit(&datastruct);
	
	if(numThreads == 0) meister(numThreads);
	else sklave();
	
	MPI_Finalize();
	return 0;
}
