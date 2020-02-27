#!/bin/bash
#
# This script submits freesurfer jobs
#
#
# Set up variables
# subject directory within BIDS structure

# change to overearching bids directory


# change to your subject list
#for subject in `cat allsubs.txt` ; do
	for subject in sub-CC110033 ; do

	echo $subject
	dataDir=/lustre/archive/q10027/CamCAN_BackUp/BIDS/
	basedir=${dataDir}/${subject}

		for num_surfs in `seq 10 1 30` ; do
				surfdir=${basedir}/surfaces/equivSurfs/${num_surfs}surfs
					if [-d ${surfdir} ]; then
					rm -R ${surfdir}
					fi
		done
done
