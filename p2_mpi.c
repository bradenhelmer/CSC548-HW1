// bthelmer Braden T Helmer
// hkambha Harish Kambhampaty
// cwkavana Colin W Kavanaugh
//
// CSC 548 Parallel Systems HW1 Problem 2
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* first grid point */
#define XI 1.0
/* last grid point */
#define XF 100.0

// Root Proc
#define ROOTPROC 0

/* function declarations */
double fn(double);
void print_y_data(int, double *, double *, double *);
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
  slice_size = ceil(((double)NGRID / numproc));
  x_slice = (double *)malloc(sizeof(double) * slice_size);
  // Allocate room for the boundaries
  y_slice = (double *)malloc(sizeof(double) * (slice_size + 2));

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

  // MESSAGE PASSING BOUNDARY VALUES
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Blocking Sends
  if (!BLOCKING) {
    if (rank < (numproc - 1)) {
      // Forward send
      MPI_Send(y_slice + slice_size, 1, MPI_DOUBLE, rank + 1, 0,
               MPI_COMM_WORLD);
      // Receive from next
      MPI_Recv(y_slice + slice_size + 1, 1, MPI_DOUBLE, rank + 1, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank > ROOTPROC) {
      // Backward Send
      MPI_Send(y_slice + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
      // Receive from previous
      MPI_Recv(y_slice, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  }
  // Non-blocking sends
  else {
    if (rank < (numproc - 1)) {
      MPI_Request requests[2];
      // Forward send
      MPI_Isend(y_slice + slice_size, 1, MPI_DOUBLE, rank + 1, 0,
                MPI_COMM_WORLD, &requests[0]);
      // Receive from next
      MPI_Irecv(y_slice + slice_size + 1, 1, MPI_DOUBLE, rank + 1, 0,
                MPI_COMM_WORLD, &requests[1]);
      MPI_Waitall(2, requests, MPI_STATUS_IGNORE);
    }
    if (rank > ROOTPROC) {
      MPI_Request requests[2];
      // Backward Send
      MPI_Isend(y_slice + 1, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                &requests[0]);
      // Receive from previous
      MPI_Irecv(y_slice, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,
                &requests[1]);
      MPI_Waitall(2, requests, MPI_STATUS_IGNORE);
    }
  }

  // DERIVATIVE CALCULATIONS
  // ~~~~~~~~~~~~~~~~~~~~~~~
  double dx = x_slice[1] - x_slice[0];
  for (loop_idx = 0; loop_idx < slice_size; loop_idx++) {
    dy_slice[loop_idx] =
        (y_slice[loop_idx + 2] - y_slice[loop_idx]) / (2.0 * dx);
  }

  // GATHERING CODE HERE FOR ROOT PROC
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Declare full grid arrays to be gathered from proc 0
  double *x_vec, *y_vec, *dy_vec;

  // allocate vectors for gathering
  if (rank == ROOTPROC) {
    x_vec = (double *)malloc(sizeof(double) * (slice_size * numproc));
    y_vec = (double *)malloc(sizeof(double) * (slice_size * numproc));
    dy_vec = (double *)malloc(sizeof(double) * (slice_size * numproc));
  }

  // MPI_Gather
  if (!GATHER_T) {
    // Root proc
    if (rank == ROOTPROC) {

      // define counts of recvs
      int recv_counts[numproc];
      // define displacemets
      int displs[numproc];
      // set values of recv_counts and displs
      for (int rank_iter = 0; rank_iter < numproc; rank_iter++) {
        recv_counts[rank_iter] = slice_size;
        displs[rank_iter] = rank_iter * slice_size;
      }

      if (!BLOCKING) {
        MPI_Gatherv(x_slice, slice_size, MPI_DOUBLE, x_vec, recv_counts, displs,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gatherv(y_slice + 1, slice_size, MPI_DOUBLE, y_vec, recv_counts,
                    displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gatherv(dy_slice, slice_size, MPI_DOUBLE, dy_vec, recv_counts,
                    displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      } else {
        // Request handler
        MPI_Request requests[3];
        MPI_Igatherv(x_slice, slice_size, MPI_DOUBLE, x_vec, recv_counts,
                     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Igatherv(y_slice + 1, slice_size, MPI_DOUBLE, y_vec, recv_counts,
                     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests[1]);
        MPI_Igatherv(dy_slice, slice_size, MPI_DOUBLE, dy_vec, recv_counts,
                     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests[2]);
        MPI_Waitall(3, requests, MPI_STATUS_IGNORE);
      }

    } else {
      // Non-root proc

      if (!BLOCKING) {
        MPI_Gatherv(x_slice, slice_size, MPI_DOUBLE, NULL, NULL, NULL,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gatherv(y_slice + 1, slice_size, MPI_DOUBLE, NULL, NULL, NULL,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gatherv(dy_slice, slice_size, MPI_DOUBLE, NULL, NULL, NULL,
                    MPI_DOUBLE, 0, MPI_COMM_WORLD);
      } else {
        // Request handler
        MPI_Request requests[3];
        MPI_Igatherv(x_slice, slice_size, MPI_DOUBLE, NULL, NULL, NULL,
                     MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Igatherv(y_slice + 1, slice_size, MPI_DOUBLE, NULL, NULL, NULL,
                     MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests[1]);
        MPI_Igatherv(dy_slice, slice_size, MPI_DOUBLE, NULL, NULL, NULL,
                     MPI_DOUBLE, 0, MPI_COMM_WORLD, &requests[2]);
        MPI_Waitall(3, requests, MPI_STATUS_IGNORE);
      }
    }
  }

  // Manual Gathering start here
  else {
    // Root proc
    if (rank == ROOTPROC) {

      // Copy root slice into vec
      memcpy(x_vec, x_slice, sizeof(double) * slice_size);
      memcpy(y_vec, y_slice + 1, sizeof(double) * slice_size);
      memcpy(dy_vec, dy_slice, sizeof(double) * slice_size);

      // Blocking manual gather
      if (!BLOCKING) {
        // Loop over remaining processes
        for (loop_idx = 1; loop_idx < numproc; loop_idx++) {
          const u_int16_t offset = loop_idx * slice_size;
          MPI_Recv(x_vec + offset, slice_size, MPI_DOUBLE, loop_idx, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(y_vec + offset, slice_size, MPI_DOUBLE, loop_idx, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(dy_vec + offset, slice_size, MPI_DOUBLE, loop_idx, 0,
                   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
      }
      // Non blocking manual gather
      else {
        const u_int8_t req_count = 3 * (numproc - 1);
        MPI_Request requests[req_count];
        for (loop_idx = 1; loop_idx < numproc; loop_idx++) {
          const u_int8_t actual_index = (loop_idx - 1) * 3;
          const u_int16_t offset = loop_idx * slice_size;
          MPI_Irecv(x_vec + offset, slice_size, MPI_DOUBLE, loop_idx, 0,
                    MPI_COMM_WORLD, &requests[actual_index]);
          MPI_Irecv(y_vec + offset, slice_size, MPI_DOUBLE, loop_idx, 0,
                    MPI_COMM_WORLD, &requests[actual_index + 1]);
          MPI_Irecv(dy_vec + offset, slice_size, MPI_DOUBLE, loop_idx, 0,
                    MPI_COMM_WORLD, &requests[actual_index + 2]);
        }
        MPI_Waitall(req_count, requests, MPI_STATUS_IGNORE);
      }
    }
    // Non root proc
    else {
      // Blocking manual send
      if (!BLOCKING) {
        MPI_Send(x_slice, slice_size, MPI_DOUBLE, ROOTPROC, 0, MPI_COMM_WORLD);
        MPI_Send(y_slice + 1, slice_size, MPI_DOUBLE, ROOTPROC, 0,
                 MPI_COMM_WORLD);
        MPI_Send(dy_slice, slice_size, MPI_DOUBLE, ROOTPROC, 0, MPI_COMM_WORLD);
      }
      // Non blocking manual send
      else {
        MPI_Request requests[3];
        MPI_Isend(x_slice, slice_size, MPI_DOUBLE, ROOTPROC, 0, MPI_COMM_WORLD,
                  &requests[0]);
        MPI_Isend(y_slice + 1, slice_size, MPI_DOUBLE, ROOTPROC, 0,
                  MPI_COMM_WORLD, &requests[1]);
        MPI_Isend(dy_slice, slice_size, MPI_DOUBLE, ROOTPROC, 0, MPI_COMM_WORLD,
                  &requests[2]);
        MPI_Waitall(3, requests, MPI_STATUS_IGNORE);
      }
    }
  }

  if (rank == ROOTPROC) {
    print_y_data(NGRID, x_vec, y_vec, dy_vec);
    free(x_vec);
    free(y_vec);
    free(dy_vec);
  }

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

  // fclose(fp);
}
