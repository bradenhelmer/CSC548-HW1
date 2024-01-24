// bthelmer Braden T Helmer
//
//
//
// CSC 548 Parallel Systems HW1 Problem 2
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/* first grid point */
#define XI 1.0
/* last grid point */
#define XF 100.0

// Root Proc
#define ROOTPROC 0

/* function declarations */
double fn(double);
void print_function_data(int, double *, double *, double *);
int main(int, char **);

int main(int argc, char *argv[]) {

  // Argument Collection
  int NGRID;
  int BLOCKING;
  int GATHER_T;
  if (argc > 3) {
    NGRID = atoi(argv[1]);
    BLOCKING = atoi(argv[2]);
    GATHER_T = atoi(argv[3]);
    if (!(BLOCKING == 0 || BLOCKING == 1)) {
      printf("Blocking value must be 0 or 1\n");
      exit(0);
    }
    if (!(GATHER_T == 0 || GATHER_T == 1)) {
      printf("Gather type value must be 0 or 1\n");
      exit(0);
    }
  } else {
    printf("Please specify grid points, blocking, and gather type values.\n");
    exit(0);
  }

  // MPI initilization
  // Number of processes and rank for current process.
  int numproc, rank;

  // Init MPI
  MPI_Init(&argc, &argv);

  // Get proc count
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);

  // Get proc rank
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Slice delcarations
  int slice_size, slice_min, slice_max;
  double *x_slice, *y_slice, *dy_slice;
  // Slice size and allocations
  slice_size = NGRID / numproc;
  x_slice = (double *)malloc(sizeof(double) * slice_size);
  // Allocate room for the boundaries
  y_slice = (double *)malloc(sizeof(double) * slice_size + 2);

  // Min and max values for domain slices.
  slice_min = rank * slice_size + 1;
  slice_max = (rank + 1) * slice_size;

  // Core calculations
  // -----------------

  // Construct domain and function slice values
  int loop_idx;
  for (loop_idx = slice_min; loop_idx <= slice_max; loop_idx++) {
    // Since the slice values are different for each process, we need calcualte
    // the actual slice index relative to the rank.
    int actual_index = loop_idx - (rank * slice_size) - 1;

    // Calculate x value
    x_slice[actual_index] =
        XI + (XF - XI) * (double)(loop_idx - 1) / (double)(NGRID - 1);

    // Calculate y value, add one index to keep in slice bounds without
    // boundaries
    y_slice[actual_index + 1] = fn(x_slice[actual_index]);
  }

  // Alloc derivative slice
  dy_slice = (double *)malloc(sizeof(double) * slice_size);

  // Statically calculate left_boundary for left edge proc
  if (rank == ROOTPROC) {
    double x0 = x_slice[1] - (x_slice[2] - x_slice[1]);
    y_slice[0] = fn(x0);
  }

  // Statically calculate right_boundary for right edge proc
  if (rank == numproc - 1) {
    double xn = x_slice[slice_size - 1] +
                (x_slice[slice_size - 1] - x_slice[slice_size - 2]);
    y_slice[slice_size + 1] = fn(xn);
  }

  // Blocking Sends
  if (!BLOCKING) {
    if (rank < numproc) {
      // Forward send
      MPI_Send(y_slice + slice_size, 1, MPI_DOUBLE, rank + 1, 0,
               MPI_COMM_WORLD);
      // Recieve from next
      MPI_Recv(y_slice + slice_size + 1, 1, MPI_DOUBLE, rank + 1, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank > ROOTPROC) {
      // Backward Send
      MPI_Send(y_slice + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      // Receive from previous
      MPI_Recv(y_slice + slice_size + 1, 1, MPI_DOUBLE, rank + 1, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
  // Non-blocking sends
  else {
    MPI_Request requests[4];
    if (rank < numproc) {
      // Forward send
      MPI_Isend(y_slice + slice_size, 1, MPI_DOUBLE, rank + 1, 0,
                MPI_COMM_WORLD, &requests[0]);
      // Recieve from next
      MPI_Irecv(y_slice + slice_size + 1, 1, MPI_DOUBLE, rank + 1, 0,
                MPI_COMM_WORLD, &requests[1]);
    }
    if (rank > ROOTPROC) {
      // Backward Send
      MPI_Isend(y_slice + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                &requests[2]);
      // Receive from previous
      MPI_Irecv(y_slice + slice_size + 1, 1, MPI_DOUBLE, rank + 1, 0,
                MPI_COMM_WORLD, &requests[3]);
    }
	MPI_Waitall(4, requests, MPI_STATUS_IGNORE);
  }

  // GATHERING CODE HERE FOR ROOT PROC
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Declare full grid arrays to be gathered from proc 0
  double *x_vec, *y_vec, *dy_vec;

  // Free proc slice arrays
  free(x_slice);
  free(y_slice);
  free(dy_slice);

  // Exit MPI nicely by finalizing.
  MPI_Finalize();
  return 0;
}

// prints out the y and its dy to a file
void print_y_data(int np, double *x, double *y, double *dydx) {
  int i;

  char filename[1024];
  sprintf(filename, "fn-%d.dat", np);

  FILE *fp = fopen(filename, "w");

  for (i = 0; i < np; i++) {
    fprintf(fp, "%f %f %f\n", x[i], y[i], dydx[i]);
  }

  fclose(fp);
}
