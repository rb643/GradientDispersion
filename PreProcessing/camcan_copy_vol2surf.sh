#!/bin/bash

# change to overearching bids directory
baseDir=/lustre/archive/q10027/
vol2surfDir=${baseDir}/CamCAN_BackUp/vol2surf/
mkdir $vol2surfDir
# read your subject list
for sub in `cat CC_subjects.txt` ; do
#for sub in CC110033 ; do

  if [ -e ${baseDir}/CamCAN_BackUp/BIDS/${sub}/proc_rsfmri/surfaces/${sub}_fmri2fs_bbr_lh_c69-32k.mgh ]
  then

    echo 'copying subject ' ${sub}
	# change to conform with your folder structure
	cp ${baseDir}/CamCAN_BackUp/BIDS/${sub}/proc_rsfmri/surfaces/${sub}_fmri2fs_bbr_lh_c69-32k.mgh ${vol2surfDir}/
	cp ${baseDir}/CamCAN_BackUp/BIDS/${sub}/proc_rsfmri/surfaces/${sub}_fmri2fs_bbr_rh_c69-32k.mgh ${vol2surfDir}/
  fi

done
