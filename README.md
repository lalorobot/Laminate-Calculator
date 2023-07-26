# Laminate-Calculator
A calculator for laminated composite plate mechanics written in Fortran.

It should be noted that this calculator is dependent on the csv_file module found within FLIBS.

compMat.f90 contains all the required subroutines for the calculation of the Q, Q̅, ABD matrices as well as the laminate's equivalent moduli and the material's plane properties.

planeProps.f90 is the main laminate setup file, it determines the fibre and matrix properties as well as the number of layers, desired angles and layer thicknesses. It prints each step of the calculation to the terminal as well as to a CSV file for further manipulation.

To-Do:

1. Automate finding the ̅z factors, right now the loop requires modifications should the number of layers change.
2. Find the ABD matrix's inverse abd.
3. Separate the individual calculation steps and the final ABD matrix and properties into two separate CSV files.

How to compile:
gfortran csv_file.f90 compMat.f90 planeProps.f90 -o planeProps
