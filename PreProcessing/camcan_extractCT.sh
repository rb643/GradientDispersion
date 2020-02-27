#!/bin/bash
#
# Set up variables
module load freesurfer

# subject directory within BIDS structure
subject=$1

# change to overearching bids directory
topDir=/lustre/archive/q10021/Cam_CAN/BIDS/
baseDir=${topDir}/${subject}/

export SUBJECTS_DIR="$baseDir"/surfaces

#mris_anatomical_stats -a "$baseDir"/surfaces/${subject}/label/lh.sjh.annot -f "$baseDir"/surfaces/${subject}/stats/lh.sjh.stats -b ${subject} lh
#mris_anatomical_stats -a "$baseDir"/surfaces/${subject}/label/rh.sjh.annot -f "$baseDir"/surfaces/${subject}/stats/rh.sjh.stats -b ${subject} rh

aparcstats2table --subjects ${subject} --hemi lh --meas thickness --parc HCP.fsaverage.aparc --tablefile "$baseDir"/surfaces/${subject}/stats/lh_${subject}_thickness_HCP.csv
aparcstats2table --subjects ${subject} --hemi rh --meas thickness --parc HCP.fsaverage.aparc --tablefile "$baseDir"/surfaces/${subject}/stats/rh_${subject}_thickness_HCP.csv

cp "$baseDir"/surfaces/${subject}/stats/lh_${subject}_thickness_sjh.csv /lustre/archive/q10021/Cam_CAN/DataOut/CT_HCP_lh/lh_${subject}_thickness_HCP.csv
cp "$baseDir"/surfaces/${subject}/stats/rh_${subject}_thickness_sjh.csv /lustre/archive/q10021/Cam_CAN/DataOut/CT_HCP_rh/rh_${subject}_thickness_HCP.csv
