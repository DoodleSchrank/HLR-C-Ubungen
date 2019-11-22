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
	int i;

	// wait for nu.. err data
	for(i = 1; i < numThreads; i++)
	{
		MPI_Recv(&data, 256, MPI_CHAR, i, 0, MPI_COMM_WORLD , MPI_STATUS_IGNORE);
		printf("%s\n", data);
	}

	// push threads over the brink
	char data2 = 't';
	for(i = 1; i < numThreads; i++)
	{
		MPI_Send(&data2, 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
	}
}

void subordinat()
{
	// define some integral variables
	struct timeval tv;
	char data[256], time[256], hostname[HOST_NAME_MAX], data2;

	// get hostname and time of day -- with failsafe, yeah!
	if(gethostname(hostname, sizeof(hostname)) != 0) strcat(hostname, "Error");
	if(gettimeofday(&tv, NULL) != 0) tv.tv_usec = 0;
	
	// string magic, wuuw
	strcat(hostname, ": ");
	sprintf(time, "%ld", tv.tv_sec);
	strcat(data, hostname);
	strcat(data, time);
	
	// send nu.. err data
	MPI_Send(&data, 256, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	
	// wait for the last push towards the brink
	MPI_Recv(&data2, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	// get ID and # of Threads for later use
	int myid, numThreads;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
	
	// Das klassische Meister-Subordinaten Prinzip
	if(myid != 0) meister(numThreads);
	else subordinat();
	
	MPI_Finalize();
	
	// build string with stringmagic
	char string[256] = "";
	if (myid == 0) snprintf(string, sizeof(string), "Meister %d", myid);
	else snprintf(string, sizeof(string),  "Sklave %d", myid);
	printf("All those moments will be lost in time, like tears in rain. Time to die. (%s)\n", string);
	return 0;
}
