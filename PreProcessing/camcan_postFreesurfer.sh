#!/bin/bash

    # Set basic parameters.
    baseDirectory=$1
    mmResolution=$2
    nativeT1Regex=$3
    parcDirectory=$4
    templateDirectory=$5
    
    subject=`basename $baseDirectory`
    subjectDirectory="$baseDirectory"/surfaces/"$subject"
    conte69Directory="$baseDirectory"/surfaces/conte69/
    warpDirectory="$baseDirectory"/xfms/
    processingDirectory="$baseDirectory"/tmpProcessingPostFreesurfer/
    finalDirectory="$baseDirectory"/proc_struct
    
    # Make directories - exit if processing directory already exists (to prevent deletion of existing files at the end of this script). 
    for x in "$conte69Directory" "$processingDirectory" "$volumeDirectory"; do
        [[ ! -d "$x" ]] && mkdir -p "$x"
    done
    
    ## Deal with FreeSurfer c_ras offset 
    MatrixX=`mri_info "$subjectDirectory"/mri/brain.finalsurfs.mgz | grep "c_r" | cut -d "=" -f 5 | sed s/" "/""/g`
    MatrixY=`mri_info "$subjectDirectory"/mri/brain.finalsurfs.mgz | grep "c_a" | cut -d "=" -f 5 | sed s/" "/""/g`
    MatrixZ=`mri_info "$subjectDirectory"/mri/brain.finalsurfs.mgz | grep "c_s" | cut -d "=" -f 5 | sed s/" "/""/g`
    echo "1 0 0 ""$MatrixX" > "$subjectDirectory"/mri/c_ras.mat
    echo "0 1 0 ""$MatrixY" >> "$subjectDirectory"/mri/c_ras.mat
    echo "0 0 1 ""$MatrixZ" >> "$subjectDirectory"/mri/c_ras.mat
    echo "0 0 0 1" >> "$subjectDirectory"/mri/c_ras.mat
    
    for hemisphere in l r; do 
        # Build the conte69-32k sphere and midthickness surface
        wb_shortcuts -freesurfer-resample-prep \
                    "$subjectDirectory"/surf/"$hemisphere"h.white \
                    "$subjectDirectory"/surf/"$hemisphere"h.pial \
                    "$subjectDirectory"/surf/"$hemisphere"h.sphere.reg \
                    "$templateDirectory"/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemisphere^^}.sphere.32k_fs_LR.surf.gii \
                    "$subjectDirectory"/surf/"$hemisphere"h.midthickness.surf.gii \
                    "$conte69Directory"/"$subject"_"$hemisphere"h_midthickness_32k_fs_LR.surf.gii \
                    "$conte69Directory"/"$subject"_"$hemisphere"h_sphereReg.surf.gii
        # Resample white and pial surfaces to conte69-32k
        for surface in pial white; do
            mris_convert "$subjectDirectory"/surf/"$hemisphere"h."$surface" "$conte69Directory"/"$hemisphere"h."$surface".surf.gii
            wb_command -surface-resample \
                    "$conte69Directory"/"$hemisphere"h."$surface".surf.gii \
                    "$conte69Directory"/"$subject"_"$hemisphere"h_sphereReg.surf.gii \
                    "$templateDirectory"/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemisphere^^}.sphere.32k_fs_LR.surf.gii \
                    BARYCENTRIC \
                    "$conte69Directory"/"$subject"_"$hemisphere"h_"$surface"_32k_fs_LR.surf.gii
        done
        for surface in pial white midthickness; do
            # Apply affine transformation to pial, white, and midthickness surfaces to bring them in line with native scans. 
            wb_command -surface-apply-affine "$conte69Directory"/"$subject"_"$hemisphere"h_"$surface"_32k_fs_LR.surf.gii \
                                                    "$subjectDirectory"/mri/c_ras.mat \
                                                    "$conte69Directory"/"$subject"_"$hemisphere"h_"$surface"_32k_fs_LR_native.surf.gii     
            # Warp surfaces to MNI space 
            wb_command -surface-apply-warpfield  "$conte69Directory"/"$subject"_"$hemisphere"h_"$surface"_32k_fs_LR_native.surf.gii \
                                                        "$warpDirectory"/"$subject"_nlmni2native_"$mmResolution"mm.nii.gz \
                                                        "$conte69Directory"/"$subject"_"$hemisphere"h_"$surface"_32k_fs_LR_nlMNI152.surf.gii \
                                                        -fnirt "$warpDirectory"/"$subject"_native2nlmni_"$mmResolution"mm.nii.gz
        done        
    done
    
    ## Compute warp of native structural to Freesurfer
    export SUBJECTS_DIR="$baseDirectory"/surfaces
    bbregister --mov "$nativeT1w" --s "$subject" --reg "$warpDirectory"/$(basename $nativeT1w .nii.gz)_t1w2fs.lta --init-fsl --t1 --o "$finalDirectory"/"$subject"_t1w_"$mmResolution"mm_fsspace.nii.gz
    
    # Clear temporary files. Add some checks to be sure there's little chance of unintentional deletes. 
    find "$processingDirectory" -mindepth 1 -maxdepth 1 -type f -regextype posix-extended -regex ".*[.nii.gz]" -delete 
    
	# Remove processing directory
	find "$processingDirectory" -mindepth 0 -maxdepth 0 -type d -empty -delete
