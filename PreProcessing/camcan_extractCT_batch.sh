#!/bin/bash
#
# Set up variables
# subject directory within BIDS structure

# change to overearching bids directory
topDir=/lustre/archive/q10021/Cam_CAN/BIDS/

# change to your subject list
for subject in `cat subjectlist.txt` ; do
#for subject in sub-CC110033 ; do
subject='sub-'${subject}

	if [ -e ${topDir}/${subject}/surfaces/${subject}/label/lh.sjh.annot ]
	then

	echo '------------------------------------ working on' ${subject} '-----------------------------------'

	 sbatch --output=/lustre/archive/q10021/Cam_CAN/logs/CTextract/${subject}_CTOu2.log --nodes=1 --ntasks=1 --cpus-per-task=1 --time=00:40:00 --mem=4000 camcan_extractCT.sh ${subject}
	else

	echo '------------------------------------ no input image for ' ${subject} '-----------------------------------'

	fi

done
