#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h> 


void meister(int numThreads)
{
	char data[256];
	for(int i = 1; i < numThreads; i++)
	{
		MPI_Recv(&data, 256, MPI_CHAR, i, 0, MPI_COMM_WORLD , MPI_STATUS_IGNORE);
		printf("%s\n", data);
	}	
}

void sklave()
{
	struct timeval tv;
	char data[256], time[256], hostname[HOST_NAME_MAX];

	if(gethostname(hostname, sizeof(hostname)) != 0) strcat(hostname, "Error");
	if(gettimeofday(&tv, NULL) != 0) tv.tv_usec = 0;
	
	strcat(hostname, ": ");
	sprintf(time, "%ld", tv.tv_sec);
	strcat(data, hostname);
	strcat(data, time);

	MPI_Send(&data, 256, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int myid, numThreads;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);

	if(myid != 0) meister(numThreads);
	else sklave(myid);
	
	MPI_Finalize();
	
	char string[256] = "";
	if (myid == 0) snprintf(string, sizeof(string), "Meister %d", myid);
	else snprintf(string, sizeof(string),  "Sklave %d", myid);
	printf("All those moments will be lost in time, like tears in rain. Time to die. (%s)\n", string);
	return 0;
}
