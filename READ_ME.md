This is still temporary!

# Modular Method Computations

This repository contains MAGMA scripts used in my master's thesis on applications of the modular method to Diophantine equations of generalized Fermat type.

The main focus is the equation

x^3 + y^3 = 7^{\alpha} z^p.

## Contents

For now the repository contains:
- 7_curve.m - the code for eliminating ellitpic curves associated to the Frey curve in case of z even via B_l coefficients
- 7_curve_2.m - same, but for c odd
- 7_curve_all-m - both combined into one
- IntegralFrobeniusMatrix.m - A code used for verifying condition 3 of Diana and Nuno's elimination theorem, taken from "Integral Tate modules and splitting of primes in torsion fields of elliptic curves" by Tommaso Giorgio Centeleghe (I will probably write my own version of it at a later point, but this one has been used for now)
- val_elim.m - a code generating prime numbers p for which there can exist no solutions to the above equation for \alpha not square mod p; based on the criterium for local points introduced by Diana and Nuno. 

## Software

All scripts are written in MAGMA.

Note: GitHub may incorrectly classify `.m` files as MATLAB files. They are not MATLAB scripts.

## Reproducibility

The computations included here were used during the development of the thesis and are provided for reproducibility and reference purposes.

## Author

Monika Zakrzewska
