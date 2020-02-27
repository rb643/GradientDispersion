#!/bin/bash
#
# This script submits freesurfer jobs
#
#
# Set up variables
# subject directory within BIDS structure

# change to overearching bids directory


# change to your subject list
for subject in `cat holdout.txt` ; do
#for subject in sub-CC722891 ; do

  baseDir=/lustre/archive/q10027/CamCAN_BackUp/BIDS/${subject}

  if [ -e ${baseDir}/anat/${subject}_T2w.nii.gz ]
	then

	echo '------------------------------------ working on' ${subject} '-----------------------------------'

	 sbatch --output=/lustre/archive/q10027/CamCAN_BackUp/logs/freesurf/${subject}_freesurf.log --nodes=1 --ntasks=1 --cpus-per-task=1 --time=20:00:00 --mem=16000 camcan_run_freesurfer_sub.sh ${baseDir} ${subject}
	else

	echo '------------------------------------ no input T2 image for ' ${subject} '-----------------------------------'

	fi

done
