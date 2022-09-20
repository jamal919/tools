#!/bin/bash
# sed -i~ -e 's/\r//g' $0;diff -q $0 $0~ || { echo "repaired $0 . Running again.";exec bash $0 $@; } # auto repair windows madness

set -x # echo command lines

if [ \! "$1" ];then
  rm L2-XE-atlas-detidor.nc
  true>test.log
  bash $0 - |& tee test.log
  tar czf test.tgz test.log predictions.dat
  exit
  fi

waveList="Z0 Q1 O1 P1 K1 N2 M2 S2 K2 L2 M4 MS4"
# 2881 frames from 2008/01/01 00:00:00 to 2008/04/30 00:00:00
simu=champs_Meteo.nc
topo=$simu

printf "%d\n%g %g\n%g %g\n" 2 -2 49 -5 45 >control.dat

if [ -e L2-XE-atlas-detidor.nc ];then
  rm L2-XE-atlas.nc
else
  # HARMONIC ANALYSIS AND DETIDING
  comodo-detidor $simu --only-atlases --variable-mask -v XE -c control.dat -d $waveList &
  comodo-control $simu -v XE -c control.dat
  wait
  comodo-detidor $simu --only-atlases -1 --variable-mask -v XE -d $waveList
  rm *-XE-*.txt
  comodo-detidor control.nc -v XE -c control.dat -d $waveList
  paste -d \  series-XE-*.txt >series.dat
  spectrum-test
  bash -c 'rm [1-9] [1-9]*[0-9]'
  paste constants-XE-*.txt
  comodo-detidor $simu --only-atlases -v U -d $waveList
  comodo-detidor $simu --only-atlases -v V -d $waveList
  mv L2-XE-atlas.nc L2-XE-atlas-detidor.nc
  fi

# ADMITTANCE METHOD : check list of wave names ?
comodo-admittance -a WAVE-XE-atlas.nc -v XE_a XE_G N2 M2 L2 K2
# compare atlases
# ncview L2-XE-atlas-detidor.nc L2-XE-atlas.nc

# INTERPOLATION AND PREDICTION
predictor -p control.dat -a WAVE-XE-atlas.nc -v XE_a XE_G -s 01/01/2008 -f 30/04/2008 -w $waveList
interp-test
head predictions.dat

# ENERGY BUDGET
comodo-energy $topo H0 WAVE-VAR-atlas.nc XE U V $waveList

true took ${SECONDS}s
which pocvip && exit

size="700 800"
type=png
pocvip -x \
  -f K1-XE-atlas.nc -v XE_a -s x -l '0 52.5 -5 30' -n K1-atlas.$type -p "$size" \
  -f M2-XE-atlas.nc -v XE_a -s x -n M2-atlas.$type -p "$size" \
  -f M2-energy.nc -v uflux -s x -v vflux -s y -n M2-flux.$type -p "$size" \
    -v pressureWork -s x -n M2-pressureWork.$type -p "$size" \
    -v dissipation -s x -n M2-dissipation.$type -p "$size" \
    -v cinE -s x -n M2-cinE.$type -p "$size" \
    -v divFlux -s x -n M2-divFlux.$type -p "$size" \
      -l '0 49 -5 180' -n M2-divFlux-l_49_-5_180.$type -p "$size" \
    -v h -s x -f /home/softs/pocvip/data/palette.txt -s "p0 180" -n bathy.$type -p "$size" \
  -q
