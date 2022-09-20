#! /bin/bash
#SBATCH --output="detidor-%j"

# ------------------------------------------------------------
# PARAMETERS

# Input and output tides series
DATA_PATH=/gsa/fles/gloss/hourly-20130416
IN=$DATA_PATH/origin
OUT=$DATA_PATH

# WAVES and ATLAS selection
FES=/gsa/fles/FES2004/netcdf/WAVE.tide.nc
ONDES=/gsa/fles/ondes/ondes65.dut

# starting year
START=1980
#START=2000
#END=1999
STEP_YEAR=2

# main program
DETIDOR=/home2/fles/tugo/src/tools/objects/detidor

# ------------------------------------------------------------
# LET'S GO

# Main command with common parameters
# -N <nb_nodes>  -n <nb_tasks>
#RUN="srun -N10 -n1 $DETIDOR -spectrum $ONDES -tides $FES -step ${STEP_YEAR}y  -s $START -context ctoh "
#RUN="srun -N5 -n10  $DETIDOR -spectrum $ONDES -tides $FES -step ${STEP_YEAR}y  -s $START -e $END -context ctoh "
RUN="srun -N5 -n10  $DETIDOR -spectrum $ONDES -tides $FES -step ${STEP_YEAR}y  -s $START -context ctoh "

echo "BEGIN"
hostname

# dispatch the list of stations into several runs

# # test
# for j in 4
#   do
#   echo $RUN $IN/h$j*.dat
#   $RUN $IN/h$j*.dat &
#   sleep 1
# done

# short series : all at once
for j in 4 6 7 8 9
  do
  $RUN $IN/h$j*.dat &
  sleep 1
done

# medium series : split in 2 parts
for j in 1 2 3 5
  do
  $RUN $IN/h$j[0-4]*.dat &
  sleep 1
  
  $RUN $IN/h$j[5-9]*.dat &
  sleep 1
done

# the serie h0* is the longest : split into 5 runs
$RUN $IN/h00*.dat  &
sleep 1

$RUN $IN/h0[1-3]*.dat  &
sleep 1

$RUN $IN/h04*.dat  &
sleep 1

$RUN $IN/h05*.dat  &
sleep 1

$RUN $IN/h0[6-9]*.dat  &
sleep 1



echo "waiting"
wait
echo "END"

#sbatch /home2/fles/tugo/scripts/run_detidor.sh
#squeue


# fles@nuwa:/gsa/fles/gloss/hourly-20120625/analysis.ondes65.FES2004.01y.save>  for i in {0..9}; do echo h$i'*.stat'; ls -lt h$i*.stat|wc;done
# h0*.stat
#    1489   11912  148900
# h1*.stat
#     881    7048   88100
# h2*.stat
#     823    6584   82300
# h3*.stat
#     891    7128   89100
# h4*.stat
#      78     624    7800
# h5*.stat
#     492    3936   49200
# h6*.stat
#     101     808   10100
# h7*.stat
#     275    2200   27500
# h8*.stat
#     384    3072   38400
# h9*.stat
#      43     344    4300
# 


# fles@nuwa:/gsa/fles/gloss/hourly-20130416/analysis.ondes65.FES2004.01y> for i in {0..9}; do echo h$i'*.stat'; ls -lt h$i*.stat|wc;done
# h0*.stat
#    1544   12352  154400
# h1*.stat
#     945    7560   94500
# h2*.stat
#     855    6840   85500
# h3*.stat
#     925    7400   92500
# h4*.stat
#      84     672    8400
# h5*.stat
#     510    4080   51000
# h6*.stat
#     105     840   10500
# h7*.stat
#     318    2544   31800
# h8*.stat
#     405    3240   40500
# h9*.stat
#      54     432    5400

# for i in {0..9}; do echo h$i'*.stat'; ls -lt h${i}????????????2010*.stat|wc;done

# for i in {0..9}; do echo h${i} 2000*.stat; ls -lt h${i}????????????2000*.stat|wc;done


# for i in {2000..2013}; do echo year ${i}; ls -lt h?????????????${i}*.stat|wc;done
