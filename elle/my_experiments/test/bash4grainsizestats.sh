#!/bin/bash
#
# By Florian Steinbach, florian.steinbach@uni-tuebingen.de
#

DEBUG=
#DEBUG=echo

rm *~

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                   INPUT                                   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Using the WritePhaseGrainSizeStatistics special function of FS_statistics for phase 1 (ice).
# Output is the grain size statistics for this phase

FILEROOT=./results/egrip_manual_tidy # Name until "_step" and without the 3 digit step number, should be the .elle file
OUTROOT=$FILEROOT"_grainsizes_ice" # output filename is defined in code

STARTSTEP=1 # If you want to start with file "_step001", type 1
STEPS=10 # Maximum number of steps
INCREMENT=1

DATAFORPHASE=1 # Type 1 for ice, 2 for bubbles or whatever 2nd phase

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                   BEGIN CALCULATIONS                      #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


echo "Extracting grain size statistics..."

for (( i=$STARTSTEP; i<=$STEPS; i=$(($i+$INCREMENT)) ))
do
    step=$(printf "%03d" $(($i)))     
    echo "Step "$step" of "$STEPS
    
    $DEBUG FS_statistics -i $FILEROOT"_step"$step.elle -u 2 $DATAFORPHASE -n 

done

$DEBUG mv GrainSizesIce.txt $OUTROOT.txt
$DEBUG mv GrainSizeStats_Info.txt $OUTROOT"_INFO.txt"
