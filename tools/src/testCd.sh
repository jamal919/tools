#!/bin/bash

################################################################################
# EDIT BELOW

# Base command
baseCmd="mpirun -np 8 tugo-mpi NEA-spectral.intg"

# Output directory base name
odr="NEA-spectral-base25ML"

# User paramer name
n="Cd"

# TUGOm parameter section
ts="dissipation"

# TUGOm parameter name
tn="quadratic_friction_coeff"

# validation data base
g="/data1/DAC/DAC2017/validation/2017-09-04-copernicus/NWS/h.mgr"

# EDIT ABOVE

################################################################################
# DO NOT EDIT BELOW

x="$1"

d=`printf "$odr-$n%g" $x`
v="$d/validate.out"

echo -n "Checking presence of $v ... "

if [ \! -e "$v" ];then
  echo "does not exist: computing."
  #echo testing;exit 8
  
  delta="${n}_$x-delta.intg"
  printf "#%s\n    %27s = %49g //\n##\n" $ts $tn $x >$delta
  
  $baseCmd -delta $delta --output-path $d
  
  tides-validate \
    -a $d/WAVE.spectral.nc -v a_eta_LGP2 G_eta_LGP2 -unstructured LGP2 \
    -g $g -o $d/validate M2
else
  echo "exists: skipping computation."
  fi

sed -re 's/^ *M2 +(-?[0-9\.]+ +){6}([0-9\.]+) .*/\2/;t;d' $v
