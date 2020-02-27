#!/bin/bash
#
# Set up variables
# subject directory within BIDS structure

# change to overearching bids directory
topDir=/lustre/archive/q10027/CamCAN_BackUp/BIDS/

# change to your subject list
for subject in `cat CC_subjects.txt` ; do
#for subject in sub-CC110033 ; do

	if [ -e ${topDir}/${subject}/proc_rsfmri/volumetric/ICA_AROMA_task-Rest_bold/denoised_func_data_nonaggr.nii.gz ]
	then

	echo '------------------------------------ working on' ${subject} '-----------------------------------'

	 sbatch --output=/lustre/archive/q10027/CamCAN_BackUp/logs/subcortical/${subject}_subcortical.log --nodes=1 --ntasks=1 --cpus-per-task=1 --time=00:40:00 --mem=4000 camcan_subcortical_timeseries.sh ${subject}
	else

	echo '------------------------------------ no input image for ' ${subject} '-----------------------------------'

	fi

done
