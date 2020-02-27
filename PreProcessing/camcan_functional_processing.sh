#!/bin/bash
# Set up variables
module load freesurfer
source /home/rb643/py360/bin/activate

# Set basic parameters.
baseDirectory=$1
subject=$2
functionalRegex=$3
mmResolution=$4
surfaceTemplateDirectory=$5
aromaDir=$6

finalDirectory="$baseDirectory"/proc_rsfmri/volumetric/
processingDirectory="$baseDirectory"/tmpProcessingRestingFunctional/
funcDirectory="$baseDirectory"/func/
warpDirectory="$baseDirectory"/xfms/
structuralDirectory="$baseDirectory"/proc_struct/
export SUBJECTS_DIR="$baseDirectory"/surfaces/
fsDirectory="$baseDirectory"/surfaces/"$subject"/
finalSurfaceDirectory="$baseDirectory"/proc_rsfmri/"$subject"/
conte69Directory="$baseDirectory"/surfaces/conte69/

# Make directories
for x in "$finalDirectory" "$warpDirectory" "$processingDirectory" "$finalSurfaceDirectory"; do
	[[ ! -d "$x" ]] && mkdir -p "$x"
done

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
			 3dresample -orient RPI -prefix "$processingDirectory"/"$tag"_reorient -inset "$processingDirectory"/"$tag"_trDrop.nii.gz
		else
			 3dresample -orient RPI -prefix "$processingDirectory"/"$tag"_reorient -inset "$rawNifti"
		fi        
		
		## Remove slices to make an even number of slices in all directions (requisite for topup).
		dimensions=`fslhd "$processingDirectory"/"$tag"_reorient.nii.gz | grep -E "^dim[1-3]" | awk '{print $NF}'`
		newDimensions=`echo $dimensions | awk '{for(i=1;i<=NF;i++){$i=$i-($i%2);print $i}}'`
		 fslroi "$processingDirectory"/"$tag"_reorient.nii.gz "$processingDirectory"/"$tag"_sliceCut.nii.gz `echo $newDimensions | sed 's/ / 0 /g' | sed 's/^/0 /'` 0 -1
		
		## framewise displacement
		 fsl_motion_outliers -i "$processingDirectory"/"$tag"_sliceCut.nii.gz -o "$finalDirectory"/"$subject"_"$tag"_spikeRegressors_FD.1D -s "$finalDirectory"/"$subject"_"$tag"_metric_FD.1D --fd
	
		#if [[ ! -f "$finalDirectory"/ICA_AROMA_"$tag"/denoised_func_data_nonaggr.nii.gz ]] ; then
		
			## Motion correction and spikes 
			 fslmaths "$processingDirectory"/"$tag"_sliceCut.nii.gz -Tmean "$processingDirectory"/"$tag"_sliceCutMean.nii.gz
			 3dvolreg -Fourier -twopass -base "$processingDirectory"/"$tag"_sliceCutMean.nii.gz -zpad 4 -prefix "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -1Dfile "$finalDirectory"/"$subject"_"$tag".1D "$processingDirectory"/"$tag"_sliceCut.nii.gz
			 fsl_motion_outliers -i "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -o "$finalDirectory"/"$subject"_"$tag"_spikeRegressors_REFRMS.1D -s "$finalDirectory"/"$subject"_"$tag"_metric_REFRMS.1D --refrms --nomoco 
			
			## Smoothing
			 fslmaths "$finalDirectory"/"$subject"_"$tag"_fmrispace.nii.gz -kernel gauss 2.1233226 -fmean "$finalDirectory"/"$subject"_"$tag"_fmrispace_smoothed.nii.gz	
				
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
				python "$aromaDir"/ICA_AROMA.py -in "$finalDirectory"/"$subject"_"$tag"_fmrispace_smoothed.nii.gz -out "$finalDirectory"/ICA_AROMA_"$tag"/ -affmat "$warpDirectory"/"$subject"_"$tag"_fmri2t1.mat -warp "$warpDirectory"/"$subject"_native2nlmni_1mm.nii.gz -mc "$finalDirectory"/"$subject"_"$tag".1D
			else
				python "$aromaDir"/ICA_AROMA.py -in "$finalDirectory"/"$subject"_"$tag"_fmrispace_smoothed.nii.gz -out "$finalDirectory"/ICA_AROMA_"$tag"/ -affmat "$warpDirectory"/"$subject"_"$tag"_fmri2t1.mat -warp "$warpDirectory"/"$subject"_native2nlmni_1mm.nii.gz -mc "$finalDirectory"/"$subject"_"$tag".1D "$finalDirectory"/ICA_AROMA_singleecho/melodic.ica -overwrite
			fi

		
	#if [[ ! -f "$finalSurfaceDirectory"/"$subject"_"$tag"_fmri2fs_bbr_"$x"_c69-32k.mgh ]] ; then
		
		## Register to Freesurfer space
		 fslmaths "$finalDirectory"/ICA_AROMA_"$tag"/denoised_func_data_nonaggr.nii.gz -Tmean "$processingDirectory"/"$subject"_"$tag"_denoised_mean.nii.gz
		 bbregister --s "$subject" --mov "$processingDirectory"/"$subject"_"$tag"_denoised_mean.nii.gz --reg "$warpDirectory"/"$subject"_"$tag"_fmri2fs_bbr.lta --init-fsl --bold 
		
		# Register to surface
		for x in lh rh; do 
		[[ $x == lh ]] && hemisphere=l || hemisphere=r
		 mri_vol2surf \
					--mov "$finalDirectory"/ICA_AROMA_"$tag"/denoised_func_data_nonaggr.nii.gz \
					--reg "$warpDirectory"/"$subject"_"$tag"_fmri2fs_bbr.lta \
					--projfrac-avg 0.2 0.8 0.1 \
					--trgsubject "$subject" \
					--interp trilinear \
					--hemi "$x" \
					--out "$finalSurfaceDirectory"/"$subject"_"$tag"_fmri2fs_bbr_"$x".mgh

			 mri_convert "$finalSurfaceDirectory"/"$subject"_"$tag"_fmri2fs_bbr_"$x".mgh "$processingDirectory"/"$subject"_"$tag"_fmri2fs_bbr_"$x".func.gii
			 wb_command -metric-resample \
				"$processingDirectory"/"$subject"_"$tag"_fmri2fs_bbr_"$x".func.gii \
				"$conte69Directory"/"$subject"_"$hemisphere"h_sphereReg.surf.gii \
				"$surfaceTemplateDirectory"/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemisphere^^}.sphere.32k_fs_LR.surf.gii \
				ADAP_BARY_AREA \
				"$processingDirectory"/"$subject"_"$tag"_fmri2fs_bbr_"$x"_c69-32k.func.gii \
				-area-surfs \
				"$fsDirectory"/surf/"$hemisphere"h.midthickness.surf.gii \
				"$conte69Directory"/"$subject"_"$hemisphere"h_midthickness_32k_fs_LR.surf.gii
			 mri_convert "$processingDirectory"/"$subject"_"$tag"_fmri2fs_bbr_"$x"_c69-32k.func.gii "$finalSurfaceDirectory"/"$subject"_"$tag"_fmri2fs_bbr_"$x"_c69-32k.mgh
		done
		
	#fi
	
done

#matlab17b -nodisplay -nodesktop -nojvm -singleCompThread -r "postRSFMRI('"$baseDirectory"', '"$subject"'); quit"

# Clear temporary files. Add some checks to be sure there's little chance of unintentional deletes. 
find "$processingDirectory" -mindepth 1 -maxdepth 1 -type f -regextype posix-extended -regex ".*[.txt|.nii.gz|.omat|_log|.1D|.gii]" -delete 

# Remove processing directory
find "$processingDirectory" -mindepth 0 -maxdepth 0 -type d -empty -delete
