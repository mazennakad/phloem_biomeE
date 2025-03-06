#!/bin/sh

FSRCS="src/datatypes.F90 \
       src/io_mod.F90 \
       src/soil.F90 \
       src/vegetation.F90 \
       src/BiomeE.F90 \
       src/main.F90"
CPPFLAGS=''
CPPFLAGS+=" -DHydro_test"
#CPPFLAGS+=" -DGrowthOFF"
CPPFLAGS+=" -DDemographyOFF"
CPPFLAGS+=" -DPhloem_test"

#gfortran $FSRCS -o ess -I/opt/local/include -L/opt/local/lib -lnetcdff
#gfortran $FSRCS -o ess -I/Users/eweng/MACPORTS/gcc49-python3/include -L/Users/eweng/MACPORTS/gcc49-python3/lib -lnetcdff
#gfortran src/datatypes.F90 src/io_mod.F90 src/soil.F90 src/vegetation.F90 src/BiomeE.F90 src/main.F90 -DHydro_test -o ess
echo 'gfortran' $FSRCS $CPPFLAGS '-o ess'
gfortran $FSRCS $CPPFLAGS -o ess
fparameter='./para_files/parameters_BCI_hydro.nml'
echo $fparameter
cat $fparameter > ./para_files/input.nml
./ess

rm ./para_files/input.nml

#Clean up compiled files
rm ess
rm esdvm.mod
rm datatypes.mod
rm io_mod.mod
rm soil_mod.mod
rm biomee_mod.mod
