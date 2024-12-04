#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, commsize;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    size_t n = 10001;
    size_t m = 10001;
    size_t my_start = rank * (n / commsize);
    size_t my_end = rank == commsize - 1 ? ((rank + 1) * (n / commsize) - 1) : n;
    double* my_matrix = malloc((my_end - my_start) * m);
    for (size_t y = 0; y <= (my_end - my_start); y++)
    {
        for (size_t x = 0; x <= m; x++)
        {
            my_matrix[y * m + x] = (rank + 1) * x;
        }
    }

    double* others = malloc((commsize - 1) *((my_end - my_start) * m));
    

    free(my_matrix);
    free(others);
    MPI_Finalize();
    return 0;
}