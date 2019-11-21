#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>

void meister(int numThreads)
{
	MPI_Datatype datastruct;
	MPI_Status status;
	int tag;
	for(int i = 1; i < numThreads; i++)
	{
		MPI_Recv(&datastruct, 1, datastruct, i, tag , MPI_COMM_WORLD , &status);
		printf("%s:%ld\n", datastruct.hostname, datastruct.time);
	}	
}

void sklave()
{
	MPI_Datatype datastruct;
	char hostname[256];
	struct timeval tv;
	struct timezone tz;
	int tag;
	
	if(gethostname(hostname, sizeof(hostname)) != 0) hostname = "Error";
	if(gettimeofday(&tv, &tz) != 0) tv.tv_usec = 0;
	
	struct data
	{
		char *host = hostname;
		long int time = tv.tv_usec;
	} data;
	MPI_Send(&data, 1, datastruct, 0, tag, MPI_COMM_WORLD);
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int myid, numThreads;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);

	// define MPI_struct
	static MPI_Datatype datastruct;
	int blocklength[] = {256, 1};
	long int displacements[] = {0, sizeof(char[256]};
	MPI_Datatype types[] = {MPI_CHAR, MPI_LONG};

	MPI_Type_struct(2, blocklength, displacements, types, &datastruct);
	MPI_Type_commit(&datastruct);
	

	if(myid == 0) meister(numThreads);
	else sklave();
	
	MPI_Finalize();
	printf("All those moments will be lost in time, like tears in rain. Time to die. (%s)", ((myid == 0) ? "Meister" : "Sklave %d", myid));
	return 0;
}
