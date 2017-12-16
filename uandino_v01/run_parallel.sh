#!/bin/bash


g++ -fopenmp -o uandino_v01.o uandino_v01.cpp `gsl-config --cflags --libs`

./uandino_v01.o #> probsTest.csv

#mpic++ -o spectra_sampling_MCMC.o spectra_sampling_MCMC.cpp `gsl-config --cflags --libs`

#mpirun -n 2 ./spectra_sampling_MCMC.o

python plotTheThing.py
