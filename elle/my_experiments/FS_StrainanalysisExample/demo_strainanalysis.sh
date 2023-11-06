#!/bin/bash
# Example for using the strain-analysis module
# by Florian Steinbach
# florian.steinbach@uni-tuebingen.de
# or:
# FlorianSteinbach@gmx.de
#
# Version 170209_1317

DEUBG= # comment for debugging
#DEBUG=echo # uncomment for debugging

rm *~

####### USER INPUT #######

INELLEFILEROOT=mysimulation # "root" of the input elle file, without "_stepXXX.elle"
                            # ATTENTION: It must be the one from the last step that should be processed
STARTSTEP=1 # 1st step to include in finite strain tensor and strain analysis
ENDSTEP=25  # Corresponingly, the last step: Make sure you use the ENDSTEP-th Elle file, too

###########################

ellefilestep=$(printf "%03d" $ENDSTEP )



# # # Post-process the data # # #

$DEBUG FS_postprocess_strainanalysis -i $INELLEFILEROOT"_step"$ellefilestep.elle -u $STARTSTEP $ENDSTEP 1 -n

# --> Actually, now you could delete the incremental position gradient tensor files of all the steps
postprocessedfile=$INELLEFILEROOT"_postprocess_strainanalysis_steps"$STARTSTEP"-"$ENDSTEP".elle"
$DEBUG mv FS_postprocess_strainanalysis.elle $postprocessedfile
$DEBUG mv FinitePosGradTensor.txt $INELLEFILEROOT"_FiniteTensors_steps"$STARTSTEP"-"$ENDSTEP".txt"



# # Preparing the plot-file # # #

$DEBUG FS_plot_strainanalysis -i $postprocessedfile -n
plotfile=$INELLEFILEROOT"_plot_strainanalysis_steps"$STARTSTEP"-"$ENDSTEP".elle"
$DEBUG mv FS_plot_strainanalysis.elle $plotfile
