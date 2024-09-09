#######
# compile_script.sh
#
# This script compiles the base library and the main code for the KS
# solution of the Khan and Thomas (2008) model extended to incorporate
# size-dependent taxes and subsidies.
#
# 'Alternative Methods for Solving Heterogeneous Firm Models'
# Stephen Terry (2015)
#
# This Version : 12/18/15
#######


FCFLAGS="-Ofast -fno-gcse -fopenmp"

echo "Compiling the base module"

gfortran ${FCFLAGS} -c -o base_lib.o base_lib.f90

echo "Compiling the main program"

gfortran ${FCFLAGS} -c -o kt_ks_subsidy.o -I . kt_ks_subsidy.f90

echo "Creating the executable"

gfortran ${FCFLAGS} base_lib.o kt_ks_subsidy.o -o kt_ks_subsidy.exe
