#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

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
}

void subordinat()
{
	// define some integral charrays
	int length = 64;
	char data[3 * length], timestr[length], hostname[length];
	struct timeval tv;

	// get hostname and time of day -- with failsafe, yeah!
	int retval = gethostname(hostname, sizeof(hostname));
	if (retval != 0) strcat(hostname, "Error");
	retval = gettimeofday(&tv, NULL);
	
	// string magic, wuuw
	strftime(timestr, length,  "%Y-%m-%d %T", localtime(&tv.tv_sec));
	snprintf(data, sizeof(data), "%s: %s.%ld\n", hostname, timestr, tv.tv_usec);
	
	// send nu.. err data
	MPI_Send(&data, 3 * length, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	// get ID and # of Threads for later use
	int myid, numThreads;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numThreads);
	
	// Das klassische Meister-Subordinaten Prinzip
	if(myid == 0) meister(numThreads);
	else subordinat();
	
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	
	// build string with stringmagic
	char string[256] = "";
	if (myid == 0) snprintf(string, sizeof(string), "Meister %d", myid);
	else snprintf(string, sizeof(string),  "Sklave %d", myid);
	printf("All those moments will be lost in time, like tears in rain. Time to die. (%s)\n", string);
	return 0;
}
