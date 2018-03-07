#!/bin/bash



g++ -fopenmp -o uandino_v02.o uandino_v02.cpp `gsl-config --cflags --libs`
while true; do {
    pgrep -x "uandino_v02" | while read PID; do {
        echo -1000 > /proc/$PID/oom_score_adj;
    } done;
} done;
#./uandino_v02.o #> probsTest.csv

#mpic++ -o spectra_sampling_MCMC.o spectra_sampling_MCMC.cpp `gsl-config --cflags --libs`

#mpirun -n 2 ./spectra_sampling_MCMC.o

python plotTheThing.py
rm uandino_v02.o
