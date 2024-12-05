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
    size_t my_start = rank == 0 ? 1 : rank * (n / commsize);
    size_t my_end = rank == commsize - 1 ? ((rank + 1) * (n / commsize) - 1) : n;
    double* my_matrix = malloc((my_end - my_start) * m);
    for (size_t y = 0; y <= (my_end - my_start); y++)
    {
        for (size_t x = 0; x <= m; x++)
        {
            my_matrix[y * m + x] = (rank + 1) * x;
        }
    }
    double* others = malloc((commsize - 1) * ((my_end - my_start) * m));
    double* ans_matr = malloc((commsize - 1) * (my_end - my_start));

	MPI_Scatter(&my_matrix, m, MPI_DOUBLE, &others, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (size_t y = my_start; y <= my_end; y++)
    {
		double mul;
		size_t x = y;
		while (my_matrix[y * m + x] == 0)
			x++;
		mul = my_matix[y * m + x] / others[y * m + x];
		for (;x <= m; x++)
		{
			my_matrix[y * m + x] = my_matrix[y * m + x] - others[y * m + x] * mul;
		}
		MPI_Scater(my_matrix[y], m, MPI_DOUBLE, &others, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	   

    free(my_matrix);
    free(others);
    MPI_Finalize();
    return 0;
}
