#######
# compile_script.sh
#
# This script compiles the base library and the main code for the solution
# of the Khan and Thomas (2008) model using the Reiter (2009) method. The
# version of the code compiled here includes shifts in micro volatility
# over the business cycle.
#
# 'Alternative Methods for Solving Heterogeneous Firm Models'
# Stephen Terry (2017)
#
# This Version : 01/16/17
#######

FCFLAGS="-O3 -fopenmp -m64"
PROJNAME=kt_reiter_vol
COMP=gfortran

echo "Removing old files."

rm base_lib.o ${PROJNAME}.o

echo "Compiling the base module."

${COMP} ${FCFLAGS} -c -o base_lib.o base_lib.f90

echo "Compiling the main program."

${COMP} ${FCFLAGS} -c -o ${PROJNAME}.o -I . ${PROJNAME}.f90

echo "Creating the executable."

${COMP} ${FCFLAGS} base_lib.o ${PROJNAME}.o -L. -lmacblas -o ${PROJNAME}.exe
