# UANdINO
Calculates the transition/survival probabilities of a(n) (anti)neutrino that travels through a certain density profile.

This program is based on the calculations done by Ohlsson and Snellman in

1. *Ohlsson, T., & Snellman, H. (2001). Neutrino oscillations with three flavors in matter of varying density. The European Physical Journal C, 20(3), 507–515. https://doi.org/10.1007/s100520100687*.

It also relies on GSL library.

2. *Galassi, M., Davies, J., Theiler, J., Gough, B., Jungman, G., Alken, P., … Ulerich, R. (2017). GNU Scientific Library Release 2.4. Retrieved from https://www.gnu.org/software/gsl/doc/latex/gsl-ref.pdf*
## Version

+ uandino_v01: Reproduces figures 1-5 from reference 1.
+ uandino_v02: Reproduces *all* figures of reference 1. Need A LOT (100k) steps to converge to a solution for figure 8, thus, it takes a lot of time.
+ uandino_py: Python implementation of the algorithm. Converges, with less steps, to the same result as v02 (1k). But still takes a lot of time to do so.
+ uandino_m: Matlab implementation. VERY fast but extremely sensitive to numerical errors, specially regarding the neutrino energy, it breaks fast at low energies.
