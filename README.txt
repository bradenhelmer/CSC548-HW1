CSC548-HW1 Group Problem

GROUP MEMBERS
-------------
bthelmer Braden T Helmer
hkambha Harish Kambhampaty
cwkavana Colin W Kavanaugh

BUILD AND RUNNING
-----------------
To build the parallel verison, run:

'make -f p2.Makefile'

To run the parallel version, run:

'mpirun -n<PROC> ./p2 <GRID POINTS> <BLOCK OPT> <GATHER OPT>'

To build the serial version, run:

'make -f p2.Makefile'

To run the serial verison, run:

'./p2_serial <GRID_POINTS>'

TESTING ANALYSIS
----------------
Criteria -
For testing the 4 scenarios (Blocking/MPI_Gather, Blocking/ManualGather, Non-blocking/MPI_Gather, Blocking/ManualGather) we tested on 2, 4, 6, 8, and 10 cores at 1 million gridpoints each, amounting to 20 different measurments.

Findings -
Testing the serial version first on 1 million gridpoints, it clocked in at 0.077 seconds. When testing the parallel version, we used a batch script that ran each scenario with a set amount of cores. We found that 2 cores provided a tremendous speedup in each of the scenarios. However, past 4 cores the speedup stalled and decreased. We believe this to be the problem of the MPI_Gather overhead. Additionally we found that the fastest scenario in each of the core configurations was using non-blockingpoint-to-point messaging along with manual gathering, reinforcing our speculation of the MPI_Gather overhead issue.
