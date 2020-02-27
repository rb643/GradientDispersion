#!/bin/bash
# Set up variables
module load freesurfer
source /home/rb643/py360/bin/activate

# Set basic parameters.
baseDirectory=$1
subject=$2
functionalRegex=$3
mmResolution=$4
templateDirectory=$5
aromaDir=$6

anatdir="$baseDirectory"/anat
finalDirectory="$baseDirectory"/proc_rsfmri/volumetric/
processingDirectory="$baseDirectory"/tmpProcessingRestingFunctional/
funcDirectory="$baseDirectory"/func
warpDirectory="$baseDirectory"/xfms/
structuralDirectory="$baseDirectory"/proc_struct/
export SUBJECTS_DIR="$baseDirectory"/surfaces/
fsDirectory="$baseDirectory"/surfaces/"$subject"/
finalSurfaceDirectory="$baseDirectory"/proc_rsfmri/"$subject"/
conte69Directory="$baseDirectory"/surfaces/conte69/

# Make directories
for x in "$finalDirectory" "$warpDirectory" "$processingDirectory" "$finalSurfaceDirectory" "$structuralDirectory"; do
	[[ ! -d "$x" ]] && mkdir -p "$x"
done

echo "running BET"
# Create reference structural image
fsl5.0-bet "$anatdir"/"$subject"_T1w.nii.gz "$warpDirectory"/bet.nii.gz  -B -f 0.1
fsl5.0-bet "$warpDirectory"/bet.nii.gz "$structuralDirectory"/"$subject"_t1w_1mm_native_brain.nii.gz -R -f 0.3

echo "running fast"
fsl5.0-fast "$structuralDirectory"/"$subject"_t1w_1mm_native_brain.nii.gz

echo "Register to MNI space and store warps." 
            fsl5.0-flirt  -ref "$templateDirectory"/MNI152_T1_1mm_brain.nii.gz -in "$structuralDirectory"/"$subject"_t1w_1mm_native_brain.nii.gz -out "$structuralDirectory"/"$subject"_t1w_1mm_lMNI152_brain.nii.gz -omat "$warpDirectory"/"$subject"_native2lmni_1mm.omat -cost mutualinfo -searchcost mutualinfo -dof 12 -interp trilinear
            fsl5.0-fnirt --ref="$templateDirectory"/MNI152_T1_1mm_brain.nii.gz --in="$structuralDirectory"/"$subject"_t1w_1mm_lMNI152_brain.nii.gz --fout="$structuralDirectory"/"$subject"_lmni2nlmni_1mm.nii.gz --interp=linear --refmask="$templateDirectory"/MNI152_T1_1mm_brain_mask.nii.gz
            
            # Create warps
           fsl5.0-convertwarp -m "$warpDirectory"/"$subject"_native2lmni_1mm.omat -w "$structuralDirectory"/"$subject"_lmni2nlmni_1mm.nii.gz -r "$templateDirectory"/MNI152_T1_1mm_brain.nii.gz -o "$structuralDirectory"/"$subject"_native2nlmni_1mm.nii.gz
           fsl5.0-invwarp -w "$structuralDirectory"/"$subject"_native2nlmni_1mm.nii.gz -o "$structuralDirectory"/"$subject"_nlmni2native_1mm.nii.gz -r "$structuralDirectory"/"$subject"_t1w_1mm_native_brain.nii.gz


# Get scans to process
restScans=`find "$funcDirectory" -maxdepth 1 -regextype posix-extended -regex "$functionalRegex".nii.gz`
echo $restScans

# Loop over all scans for everything before motion correction across scans.
for rawNifti in ${restScans}; do

	# Get basic parameters
	tag=`echo "$rawNifti" | sed "s/.*${subject}_//" | sed 's/.nii.gz//'`

	if [[ ! -f "$finalDirectory"/"$subject"_"$tag"_fmrispace_smoothed.nii.gz ]] ; then

		### Drop first five TRs and reorient.
		if [[ `echo "$rawNifti" | grep -E "$functionalRegex"` ]]; then
			nifti_tool -cbl -prefix "$processingDirectory"/"$tag"_trDrop.nii.gz -infiles "$rawNifti"'[5..$]'
			 3dresample -orient RPI -prefix "$processingDirectory"/"$tag"_reorient.nii.gz -inset "$processingDirectory"/"$tag"_trDrop.nii.gz
		else
			 3dresample -orient RPI -prefix "$processingDirectory"/"$tag"_reorient -inset "$rawNifti"
		fi

		## Remove slices to make an even number of slices in all directions (requisite for topup).
		dimensions=`fslhd "$processingDirectory"/"$tag"_reorient.nii.gz | grep -E "^dim[1-3]" | awk '{print $NF}'`
		newDimensions=`echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}'`
		 fslroi "$processingDirectory"/"$tag"_reorient.nii.gz "$processingDirectory"/"$tag"_sliceCut.nii.gz `echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /'` 0 -1
		echo "computing FD"
		## framewise displacement
		 fsl_motion_outliers -i "$processingDirectory"/"$tag"_sliceCut.nii.gz -o "$finalDirectory"/"$subject"_"$tag"_spikeRegressors_FD.1D -s "$finalDirectory"/"$subject"_"$tag"_metric_FD.1D --fd

		#if [[ ! -f "$finalDirectory"/ICA_AROMA_"$tag"/denoised_func_data_nonaggr.nii.gz ]] ; then

			## Motion correction and spikes
			 fslmaths "$processingDirectory"/"$tag"_sliceCut.nii.gz -Tmean "$processingDirectory"/"$tag"_sliceCutMean.nii.gz
			 3dvolreg -Fourier -twopass -base "$processingDirectory"/"$tag"_sliceCutMean.nii.gz -zpad 4 -prefix "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -1Dfile "$finalDirectory"/"$subject"_"$tag".1D "$processingDirectory"/"$tag"_sliceCut.nii.gz
			 fsl_motion_outliers -i "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -o "$finalDirectory"/"$subject"_"$tag"_spikeRegressors_REFRMS.1D -s "$finalDirectory"/"$subject"_"$tag"_metric_REFRMS.1D --refrms --nomoco
		echo "smoothing"
			## Smoothing
			 fslmaths "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -kernel gauss 2.1233226 -fmean "$finalDirectory"/"$subject"_"$tag"_fmrispace_smoothed.nii.gz
		
		echo "running flirt"
			## Convert structural pves to functional space
			 flirt -in "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -ref "$structuralDirectory"/"$subject"_t1w_"$mmResolution"mm_native_brain.nii.gz -omat "$warpDirectory"/"$subject"_"$tag"_fmri2t1.mat
			 convert_xfm -inverse "$warpDirectory"/"$subject"_"$tag"_fmri2t1.mat -omat "$warpDirectory"/"$subject"_"$tag"_t12fmri.mat
			for idx in 0 1 2; do
				[[ $idx == 0 ]] && tissue=CSF
				[[ $idx == 1 ]] && tissue=GM
				[[ $idx == 2 ]] && tissue=WM
				tissuemap=$(find "$structuralDirectory" -name "*native_brain_pve_${idx}.nii.gz")
				 flirt -ref "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -in "$tissuemap" -applyxfm -init "$warpDirectory"/"$subject"_"$tag"_t12fmri.mat -out "$processingDirectory"/"$subject"_"$tag"_"$tissue".nii.gz
				 fslmeants -i "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -o "$finalDirectory"/"$subject"_"$tag"_"$tissue".txt -m "$processingDirectory"/"$subject"_"$tag"_"$tissue".nii.gz -w
			done
			 fslmaths "$processingDirectory"/"$subject"_"$tag"_WM.nii.gz -add  "$processingDirectory"/"$subject"_"$tag"_GM.nii.gz -add  "$processingDirectory"/"$subject"_"$tag"_CSF.nii.gz "$processingDirectory"/"$subject"_"$tag"_WB.nii.gz
			 fslmeants -i "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -o "$finalDirectory"/"$subject"_"$tag"_global.txt -m "$processingDirectory"/"$subject"_"$tag"_WB.nii.gz -w
	fi

			## ica aroma
			if [[ ! -d "$finalDirectory"/ICA_AROMA_"$tag"/ ]] ; then
				python "$aromaDir"/ICA_AROMA.py -in "$finalDirectory"/"$subject"_"$tag"_fmrispace_smoothed.nii.gz -out "$finalDirectory"/ICA_AROMA_"$tag"/ -affmat "$warpDirectory"/"$subject"_"$tag"_fmri2t1.mat -warp "$structuralDirectory"/"$subject"_native2nlmni_1mm.nii.gz -mc "$finalDirectory"/"$subject"_"$tag".1D
			else
				python "$aromaDir"/ICA_AROMA.py -in "$finalDirectory"/"$subject"_"$tag"_fmrispace_smoothed.nii.gz -out "$finalDirectory"/ICA_AROMA_"$tag"/ -affmat "$warpDirectory"/"$subject"_"$tag"_fmri2t1.mat -warp "$warpDirectory"/"$subject"_native2nlmni_1mm.nii.gz -mc "$finalDirectory"/"$subject"_"$tag".1D "$finalDirectory"/ICA_AROMA_singleecho/melodic.ica -overwrite
			fi



done

#matlab17b -nodisplay -nodesktop -nojvm -singleCompThread -r "postRSFMRI('"$baseDirectory"', '"$subject"'); quit"

# Clear temporary files. Add some checks to be sure there's little chance of unintentional deletes.
#find "$processingDirectory" -mindepth 1 -maxdepth 1 -type f -regextype posix-extended -regex ".*[.txt|.nii.gz|.omat|_log|.1D|.gii]" -delete

# Remove processing directory
#find "$processingDirectory" -mindepth 0 -maxdepth 0 -type d -empty -delete
