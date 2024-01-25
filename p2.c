// bthelmer Braden T Helmer
// cwkavana Colin W Kavanaugh
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

  //MPI initilization requierments
  /* process information */
  int   numproc, rank, len;
  /* current process hostname */
  char  hostname[MPI_MAX_PROCESSOR_NAME];
  /* initialize MPI */
  MPI_Init(&argc, &argv);
  /* get the number of procs in the comm */
  MPI_Comm_size(MPI_COMM_WORLD, &numproc);
  /* get my rank in the comm */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /* get some information about the host I'm running on */
  MPI_Get_processor_name(hostname, &len);

  

  // loop index
  int i;

  // domain array and step size
  double *xc = (double *)malloc(sizeof(double) * (NGRID + 2));
  double dx;

  // function array and derivative
  // the size will be dependent on the
  // number of processors used
  // to the program
  double *yc, *dyc;

  //"real" grid indices
  int imin, imax;

  imin = 1;
  imax = NGRID;

  // construct grid
  for (i = 1; i <= NGRID; i++) {
    xc[i] = XI + (XF - XI) * (double)(i - 1) / (double)(NGRID - 1);
  }
  // step size and boundary points
  dx = xc[2] - xc[1];
  xc[0] = xc[1] - dx;
  xc[NGRID + 1] = xc[NGRID] + dx;

  // allocate function arrays
  yc = (double *)malloc((NGRID + 2) * sizeof(double));
  dyc = (double *)malloc((NGRID + 2) * sizeof(double));

  // define the function
  for (i = imin; i <= imax; i++) {
    yc[i] = fn(xc[i]);
  }

  // set boundary values
  yc[imin - 1] = 0.0;
  yc[imax + 1] = 0.0;

  // NB: boundary values of the whole domain
  // should be set
  yc[0] = fn(xc[0]);
  yc[imax + 1] = fn(xc[NGRID + 1]);

  // compute the derivative using first-order finite differencing
  //
  //   d           f(x + h) - f(x - h)
  //  ---- f(x) ~ --------------------
  //   dx                 2 * dx
  //
  for (i = imin; i <= imax; i++) {
    dyc[i] = (yc[i + 1] - yc[i - 1]) / (2.0 * dx);
  }

  print_function_data(NGRID, &xc[1], &yc[1], &dyc[1]);

  // free allocated memory
  free(yc);
  free(dyc);

  MPI_Finalize();
  
  return 0;
}

<<<<<<< Updated upstream
// prints out the function and its derivative to a file
void print_function_data(int np, double *x, double *y, double *dydx) {
=======
// prints out the y and its dy to a file
void print_y_data(int np, double *x, double *y, double *dydx) {
  //
>>>>>>> Stashed changes
  int i;
  

  //
  char filename[1024];
  //printf("the filename is %d\n", filename[1024]);
  //char * commandsForGnuplot[] = {"set title \"%d\"", np, "plot '%d'", filename};
  
  
  //
  //printf("Printing the filename");
  sprintf(filename, "fn-%d.dat", np);
  printf("printing file name %s \n", filename);
  
  char * commandsForGnuplot[] = {"set title \"p2\"", "plot 'fn-100.dat'"};
  
  printf("printing commands for gnuplot: %s\n", commandsForGnuplot[1]);

  //
  FILE *fp = fopen(filename, "w");
  
  
  //Opens a gnuplot file to be displayed
  FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");


  for (i = 0; i < np; i++) {
    fprintf(fp, "%f %f %f\n", x[i], y[i], dydx[i]);
  }

<<<<<<< Updated upstream
  fclose(fp);
}
=======

  //Runs the gnuplot commands and displays the grapph
  for (i=0; i < 2; i++) {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
  }
  // fclose(fp); 
  
}
>>>>>>> Stashed changes
