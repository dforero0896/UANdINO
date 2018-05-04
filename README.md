# UANdINO
Calculates the transition/survival probabilities of a(n) (anti)neutrino that travels through a certain density profile.

This program is based on the calculations done by Ohlsson and Snellman in

1. *Ohlsson, T., & Snellman, H. (2001). Neutrino oscillations with three flavors in matter of varying density. The European Physical Journal C, 20(3), 507–515. https://doi.org/10.1007/s100520100687*.

It also relies on GSL library.

2. *Galassi, M., Davies, J., Theiler, J., Gough, B., Jungman, G., Alken, P., … Ulerich, R. (2017). GNU Scientific Library Release 2.4. Retrieved from https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf*

Solar neutrino oscillations reference is

3. *Kuo, T. ., & Pantaleone, J. (1989). Neutrino oscillations in matter. Reviews of Modern Physics, 61(October), 937–979. https://doi.org/https://doi.org/10.1103/RevModPhys.61.937*
## Version

+ **uandino_v01**: Reproduces figures 1-5 from reference 1.
+ **uandino_py** and **uandino_mat** have no *skipping* scheme implemented and are located in [**legacy**](https://github.com/dforero0896/UANdINO/tree/master/legacy) folder.
+ **uandino_v02**: Reproduces *all* figures of reference 1.
### uandino_v03
Optimized through an scheme that I called *skipping*, explained in my Thesis' document, [**here**](https://github.com/dforero0896/Physics_Monograph/blob/master/physics/document/Phys_Thesis_Document.pdf). It still shows problems at low energies, nevertheless, the working threshold for the program has been moved back, significantly, to (anti)neutrino energies of about 10^3 eV. Nonetheless, the software has been tested with solar neutrinos, details on this are, again, [here](https://github.com/dforero0896/Physics_Monograph/blob/master/physics/document/Phys_Thesis_Document.pdf).

#### Running uandino:
Files `uandino.h` and `uandino.cpp` are the header and implementation files for UANdINO, therefore, they should both be located in your working directory and imported in your C++ code when intended to be used.

You may want to edit the implementation file and declare a suitable `path_resolution` variable. This defines how big the spatial step will be, the smaller, the better, if you want the program to work at low energies. For high energies, a bigger step, maybe tenths of km, should be ok. Keep in mind that if the spatial step is very small, the program will take a very long time. For small steps, you should provide some sort of `leap` greater than one, this will cause the program to skip the calculations of the operator when the density remains *almost* constant, which will be true in the small-step scale. I you do not wish to use the skipping scheme because you are using high energies, please provide `leap=1` to the function `calculateProbabilities` on calling.

Besides the `leap`, function ```void calculateProbabilities(vector<float> path, int N, int Steps, int leap, float E_min, float E_max)``` asks for a vector of floats that defines the path of length `Steps/leap`. `N` defines the number of energies to be evaluated, `Steps` the number of spatial steps to be made, something like  `total_length/path_resolution`. `E_min,max` are the energy boundaries to be used, energies are linspaced between these.
