#!/bin/bash
#
# This script run freesurfer 
#
#
# Set up variables
module load freesurfer

# subject directory within BIDS structure
basedir=$1
subject=$2
surfdir=${basedir}/surfaces/
tmpDir=/lustre/scratch/wbic-beta/rb643/temp

# set up and make necessary subfolders
for x in "$surfdir" "$tmpDir"; do
	[[ ! -d "$x" ]] && mkdir -p "$x"
done

export SUBJECTS_DIR=${surfdir}
export TMPDIR=${tmpDir}

# Run freesurufer
recon-all -subject ${subject} -i ${basedir}/anat/${subject}_T1w.nii.gz -T2 /${basedir}/anat/${subject}_T2w.nii.gz -T2pial -all -no-isrunning
