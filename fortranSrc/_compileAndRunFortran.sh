#!/bin/bash

# On Linux you need to run this first to ensure the script has execute permissions:
# chmod +x _compileAndRunFortran.sh

# This script compiles the Fortran 90 code using gfortran.
gfortran finalProjDriver.f90 finalProjMod.o -o finalProject.exe

# Run the compiled Fortran program
./finalProject.exe