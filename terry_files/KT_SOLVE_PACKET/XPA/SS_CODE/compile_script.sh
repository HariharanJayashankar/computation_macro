#######
# compile_script.sh
#
# This script compiles the base library and the main code for the XPA
# solution of the Khan and Thomas (2008) model in the steady-state.
#
# 'Alternative Methods for Solving Heterogeneous Firm Models'
# Stephen Terry (2015)
#
# This Version : 12/18/15
#######

rm base_lib.o kt_ss.o

FCFLAGS="-O3"

echo "Compiling the base module"

gfortran ${FCFLAGS} -c -o base_lib.o base_lib.f90

echo "Compiling the main program"

gfortran ${FCFLAGS} -c -o kt_ss.o -I . kt_ss.f90

echo "Creating the executable"

gfortran ${FCFLAGS} base_lib.o kt_ss.o -o kt_ss.exe