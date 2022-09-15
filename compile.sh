#compile a file in fortran
echo Give program name
read name
DFLAGS=" -Wall -fcheck=all -ffpe-trap=overflow,zero,invalid -pedantic -ggdb3"
gfortran ${DFLAGS} -o $name almo_new.f90 -llapack

