#include <mpi.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int id, num_threads;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_threads);

    int length = 128;
    char buffer[length];
    char output[length * num_threads];

    struct timeval tv;
    gettimeofday(&tv, NULL);

    char time[length];
    strftime(time, length, "%Y-%m-%d %T", localtime(&tv.tv_sec));
    char hostname[length];
    gethostname(hostname, sizeof(hostname));
    snprintf(buffer, sizeof(buffer), "%s: %s.%ld", hostname, time, tv.tv_usec);

    MPI_Gather(&buffer, length, MPI_CHAR, output, length, MPI_CHAR, 0, MPI_COMM_WORLD);

    if(id == 0 && num_threads >= 1)
    {
        for(int i = 1; i < num_threads; i++)
        {
            printf("%s\n", output + (length * i));
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Rang %d beendet jetzt!\n", id);

    MPI_Finalize();
    return 0;
}
