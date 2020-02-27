#!/bin/bash
#
# This script takes a volumetric myelin sensitive image, and evaluates
# the intensity values along precreated intracortical surfaces. Additionally,
# it will map a annotation file to the individual subject space.
#
#
# Set up variables
# subject directory within BIDS structure

# change to overearching bids directory
topDir=/lustre/archive/q10027/CamCAN_BackUp/BIDS
lhAnnot=/lustre/archive/q10027/Scripts/templates/lh.sjh.annot
rhAnnot=/lustre/archive/q10027/Scripts/templates/rh.sjh.annot
surfaceTemplateDirectory=/lustre/archive/q10027/Scripts/templates/conte69



# change to your subject list
for sub in `cat holdout.txt` ; do
#for sub in sub-CC722891 ; do
subject=${sub}

	datadir=/lustre/archive/q10027/CamCAN_BackUp/BIDS/${subject}
	echo '------------------------------------ working on' ${subject} '-----------------------------------'

	 sbatch --output=/lustre/archive/q10027/CamCAN_BackUp/logs/postFree/${subject}_func.log --nodes=1 --ntasks=1 --cpus-per-task=1 --time=04:00:00 --mem=6000 camcan_postFreesurfer.sh ${datadir} "1" ${subject} ${surfaceTemplateDirectory}


done
