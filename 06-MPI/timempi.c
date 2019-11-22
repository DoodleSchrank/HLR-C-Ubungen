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
	char data[256], timestr[256], hostname[256];

	// get hostname and time of day -- with failsafe, yeah!
	int retval = gethostname(data, sizeof(data));
	if (retval != 0) strcat(hostname, "Error");
	// get time
	time_t now;
	time(&now);
	struct tm *loctime = localtime(&now);
	
	// string magic, wuuw
	strcat(data, ": ");
	sprintf(timestr, "%s", asctime(loctime));
	strcat(data, timestr);
	
	// send nu.. err data
	MPI_Send(&data, 256, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
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
