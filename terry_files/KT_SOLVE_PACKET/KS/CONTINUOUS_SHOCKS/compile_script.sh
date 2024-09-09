#######
# compile_script.sh
#
# This script compiles the base library and the main code for the KS
# solution of the Khan and Thomas (2008) model extended to allow for
# continuous shocks.
#
# 'Alternative Methods for Solving Heterogeneous Firm Models'
# Stephen Terry (2015)
#
# This Version : 12/18/15
#######

FCFLAGS="-O3 -fopenmp -m64"
PROJNAME=kt_ks_cont
COMP=gfortran

echo "Removing old files."

rm base_lib.o ${PROJNAME}.o

echo "Compiling the base module."

${COMP} ${FCFLAGS} -c -o base_lib.o base_lib.f90

echo "Compiling the main program."

${COMP} ${FCFLAGS} -c -o ${PROJNAME}.o -I . ${PROJNAME}.f90

echo "Creating the executable."

${COMP} ${FCFLAGS} base_lib.o ${PROJNAME}.o -L. -lmacblas -lmaclapack -o ${PROJNAME}.exe