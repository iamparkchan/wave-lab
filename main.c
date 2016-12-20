#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

// Parameters
#define J 1000
#define alpha 1.
#define TIME_INTERVAL 10.
#define dt (alpha / J) // Courant Condition
#define N (TIME_INTERVAL / dt)

// Macros
#define TWO_DOT(y) (alpha * alpha * J * J * (y[i + 1] + y[i - 1] - 2. * y[i]))
#define CEIL_DIV(x, y) (((x) + (y) - 1) / (y))

void exchange(int size, int rank, double *array, int I);

int main(int argc, char *argv[])
{
    int n;
	int j;
    int i;
    double *y, *v, *a, *da, *dda;
    double *y_next, *v_next;
    int size, rank;
    int I, j0;
    MPI_File fh;
    FILE *fp;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Job Assignment
    if (rank != size - 1) {
        I = CEIL_DIV(J - 1, size) + 1;
    } else {
        I = (J - 1) - CEIL_DIV(J - 1, size) * (size - 1) + 1;
    }
    j0 = CEIL_DIV(J - 1, size) * rank;
    
    // Memory Allocation
    y      = (double*)malloc((I + 1)*sizeof(double));
    v      = (double*)malloc((I + 1)*sizeof(double));
    a      = (double*)malloc((I + 1)*sizeof(double));
    da     = (double*)malloc((I + 1)*sizeof(double));
    dda    = (double*)malloc((I + 1)*sizeof(double));
    y_next = (double*)malloc((I + 1)*sizeof(double));
    v_next = (double*)malloc((I + 1)*sizeof(double));
    
    // Record of parameters
    if (rank == 0) {
        fp = fopen("parameters.txt", "w");
        fprintf(fp, "%d %d %e", (int)N, J, TIME_INTERVAL);
        fclose(fp);
    }
    
	// Initialization of functions
	for (i = 0; i <= I; i++) {
        if ((rank == 0        && i == 0) ||
            (rank == size - 1 && i == I)) {
            y[i] = v[i] = a[i] = 0.;
        } else {
            j = j0 + i;
            if (j < J/2) {
                y[i] = pow(sin(2. * M_PI * j / J), 2.);
                v[i] = - 4. * M_PI * alpha * cos(2. * M_PI * j / J) * sin(2. * M_PI * j / J);
            } else {
                y[i] = 0.;
                v[i] = 0.;
            }
        }
	}


    // Evolution of functions
    MPI_File_open(MPI_COMM_WORLD, "output", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
    n = 0;
    do {
        exchange(size, rank, y, I);
        exchange(size, rank, v, I);
        for (i = 1; i < I; i++) {
              a[i] = TWO_DOT(y);
        }
        exchange(size, rank, a, I);
        for (i = 1; i < I; i++) {
             da[i] = TWO_DOT(v);
            dda[i] = TWO_DOT(a);
        }

        MPI_File_write_at(fh, ((J + 1) * n + j0 + ((rank == 0) ? 0 : 1)) * sizeof(double), y + ((rank == 0) ? 0 : 1), I + ((rank == 0 || rank == size - 1) ? 0 : -1), MPI_DOUBLE, MPI_STATUS_IGNORE);
        
        if (n++ >= N) break;

        for (i = 1; i < I; i++) {
            y_next[i] = y[i] + v[i] * dt +  a[i] * dt * dt / 2. +  da[i] * dt * dt * dt / 4.;
            v_next[i] = v[i] + a[i] * dt + da[i] * dt * dt / 2. + dda[i] * dt * dt * dt / 4.;
        }
        for (i = 1; i < I; i++) {
            y[i] = y_next[i];
            v[i] = v_next[i];
        }
    } while (1);
    MPI_File_close(&fh);
    
    MPI_Finalize();
    
    return 0;
}

void exchange(int size, int rank, double *array, int I)
{
    MPI_Request request;
    
    if (rank != 0       ) { MPI_Isend(array + 1    , 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &request); MPI_Request_free(&request); }
    if (rank != size - 1) { MPI_Isend(array + I - 1, 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &request); MPI_Request_free(&request); }
    if (rank != 0       ) { MPI_Recv (array + 0    , 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); }
    if (rank != size - 1) { MPI_Recv (array + I    , 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); }
}
