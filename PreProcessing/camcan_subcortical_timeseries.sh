#!/bin/bash
#
# Set up variables
module load freesurfer

subject=$1

# change to overarching bids directory
topDir=/lustre/archive/q10027/CamCAN_BackUp/BIDS/
baseDir="$topDir"/"$subject"

# set up and make necessary folders
warpDir="$baseDir"/xfms/
structDir="$baseDir"/proc_struct/
funcDir="$baseDir"/proc_rsfmri/
funcImage="$funcDir"/volumetric/ICA_AROMA_task-Rest_bold/denoised_func_data_nonaggr.nii.gz

# run first
run_first_all -b -i "$structDir"/"$subject"_t1w_1mm_native_brain.nii.gz -o "$structDir"/"$subject"_native

echo "reconstructing warping matrix"
fslmaths "$funcDir"/volumetric/"$subject"_task-Rest_bold_fmrispace.nii.gz -Tmean "$funcDir"/volumetric/"$subject"_task-Rest_bold_mean.nii.gz
fsl5.0-bet "$funcDir"/volumetric/"$subject"_task-Rest_bold_mean.nii.gz "$funcDir"/volumetric/"$subject"_task-Rest_bold_mean_bet.nii.gz -m -R -f 0.2
flirt -in "$funcDir"/volumetric/"$subject"_task-Rest_bold_mean_bet.nii.gz -ref "$structDir"/"$subject"_t1w_1mm_native_brain.nii.gz -omat "$warpDir"/"$subject"_task-Rest_bold_fmri2t1.mat 
convert_xfm -inverse "$warpDir"/"$subject"_task-Rest_bold_fmri2t1.mat -omat "$warpDir"/"$subject"_task-Rest_bold_t12fmri.mat

echo "extracting time-series"
# extract time series
for x in 10 11 12 13 16 17 18 26 49 50 51 52 53 54 58 ; do
	fslmaths "$structDir"/"$subject"_native_all_fast_firstseg.nii.gz -thr $x -uthr $x "$structDir"/"$subject"_first_$x.nii.gz
	flirt -in "$structDir"/"$subject"_first_$x.nii.gz -applyxfm -init "$warpDir"/"$subject"_task-Rest_bold_t12fmri.mat -ref "$funcImage" -out "$funcDir"/"$subject"_func_first_$x.nii.gz
	fslmeants -i $funcImage -m "$funcDir"/"$subject"_func_first_$x.nii.gz -o "$funcDir"/"$subject"_fslmeants_$x.txt
done
