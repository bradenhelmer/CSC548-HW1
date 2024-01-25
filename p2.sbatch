#!/bin/bash
# bthelmer Braden T Helmer
# hkambha Harish Kambhampaty
# cwkavana Colin W Kavanaugh

#SBATCH -J CSC548_mpi_hello     # Job name
#SBATCH -o job.%j.out           # Name of stdout output file (%j expands to jobId)
#SBATCH -N 2                    # Total number of nodes requested
#SBATCH -n 2                    # Total number of mpi tasks requested
#SBATCH -t 01:30:00             # Maximum Run time (hh:mm:ss) - 1.5 hours
#SBATCH -p rome            # Specify the broadwell queue

make -f p2.Makefile
