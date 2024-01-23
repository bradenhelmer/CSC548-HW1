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

/* function declarations */
double fn(double);
void print_function_data(int, double *, double *, double *);
int main(int, char **);

int main(int argc, char *argv[]) {

  int NGRID;
  int BLOCKING;
  int GATHER_T;
  if (argc > 3) {
    NGRID = atoi(argv[1]);
    BLOCKING = atoi(argv[2]);
    GATHER_T = atoi(argv[3]);
    if (!((BLOCKING == 0 || BLOCKING == 1) &&
          (GATHER_T == 0 || GATHER_T == 1))) {
      printf("Please specify blocking and gather values as 0 or 1\n");
      exit(0);
    }
  } else {
    printf("Please specify grid points, blocking, and gather type values.\n");
    exit(0);
  }

  // MPI initilization requierments
  /* process information */
  int numproc, rank, len;
  /* current process hostname */
  char hostname[MPI_MAX_PROCESSOR_NAME];
  /* initialize MPI */
  MPI_Init(&argc, &argv);
  /* get the number of procs in the comm */
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  /* get my rank in the comm */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /* get some information about the host I'm running on */
  MPI_Get_processor_name(hostname, &len);

  // loop index
  int loop_idx;

  // domain array and step size
  double *domain_array = (double *)malloc(sizeof(double) * (NGRID + 2));
  double step_size;

  // function array and derivative
  // the size will be dependent on the
  // number of processors used
  // to the program
  double *function_array, *derivative_array;

  //"real" grid indices
  int imin, imax;

  imin = 1;
  imax = NGRID;

  // construct grid
  for (loop_idx = 1; loop_idx <= NGRID; loop_idx++) {
    domain_array[loop_idx] =
        XI + (XF - XI) * (double)(loop_idx - 1) / (double)(NGRID - 1);
  }
  // step size and boundary points
  step_size = domain_array[2] - domain_array[1];
  domain_array[0] = domain_array[1] - step_size;
  domain_array[NGRID + 1] = domain_array[NGRID] + step_size;

  // allocate function arrays
  function_array = (double *)malloc((NGRID + 2) * sizeof(double));
  derivative_array = (double *)malloc((NGRID + 2) * sizeof(double));

  // define the function
  for (loop_idx = imin; loop_idx <= imax; loop_idx++) {
    function_array[loop_idx] = fn(domain_array[loop_idx]);
  }

  // set boundary values
  function_array[imin - 1] = 0.0;
  function_array[imax + 1] = 0.0;

  // NB: boundary values of the whole domain
  // should be set
  function_array[0] = fn(domain_array[0]);
  function_array[imax + 1] = fn(domain_array[NGRID + 1]);

  // compute the derivative using first-order finite differencing
  //
  //   d           f(x + h) - f(x - h)
  //  ---- f(x) ~ --------------------
  //   dx                 2 * dx
  //
  for (loop_idx = imin; loop_idx <= imax; loop_idx++) {
    derivative_array[loop_idx] =
        (function_array[loop_idx + 1] - function_array[loop_idx - 1]) /
        (2.0 * step_size);
  }

  print_function_data(NGRID, &domain_array[1], &function_array[1],
                      &derivative_array[1]);

  // free allocated memory
  free(function_array);
  free(derivative_array);

  return 0;
}

// prints out the function and its derivative to a file
void print_function_data(int np, double *x, double *y, double *dydx) {
  int i;

  char filename[1024];
  sprintf(filename, "fn-%d.dat", np);

  FILE *fp = fopen(filename, "w");

  for (i = 0; i < np; i++) {
    fprintf(fp, "%f %f %f\n", x[i], y[i], dydx[i]);
  }

  fclose(fp);
}
