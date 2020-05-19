#!/bin/bash
mpirun -np 16 ./jacobi 512 20000 -16 1 --mca opal_warn_on_missing_libcuda 0
