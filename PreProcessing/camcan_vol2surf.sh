#!/bin/bash
#
# Set up variables
module load freesurfer

# subject directory within BIDS structure
subject=$1

# change to overearching bids directory
topDir=/lustre/archive/q10027/CamCAN_BackUp/BIDS/
templateDir=//lustre/archive/q10027/Scripts/templates/

baseDir=${topDir}/${subject}/

# set up and make necessary subfolders
inDir="$baseDir"/proc_rsfmri/volumetric/ICA_AROMA_task-Rest_bold/
outDir="$baseDir"/proc_rsfmri/surfaces/
tmpProcessingDirectory="$baseDir"/tmpProcessingRestingFunctional
warpDir="$baseDir"/xfms
for thisDir in $warpDir $outDir $tmpProcessingDirectory; do
        [[ ! -d "$thisDir" ]] && mkdir "$thisDir"
done

export SUBJECTS_DIR="$baseDir"/surfaces
export TMPDIR=/lustre/scratch/wbic-beta/rb643/temp

for x in lh rh; do
    [[ $x == lh ]] && hemisphere=l || hemisphere=r
    mri_vol2surf \
                --mov "$inDir"denoised_func_data_nonaggr.nii.gz \
                --reg "$warpDir"/"$subject"_task-Rest_bold_fmri2fs_bbr.lta \
                --projfrac-avg 0.2 0.8 0.1 \
                --trgsubject "$subject" \
                --interp trilinear \
                --hemi "$x" \
                --out "$outDir"/"$subject"_fmri2fs_bbr_"$x".mgh

    mri_convert "$outDir"/"$subject"_fmri2fs_bbr_"$x".mgh "$tmpProcessingDirectory"/"$subject"_fmri2fs_bbr_"$x".func.gii

    wb_shortcuts -freesurfer-resample-prep \
                    "$SUBJECTS_DIR"/"$subject"/surf/"$x".white \
                    "$SUBJECTS_DIR"/"$subject"/surf/"$x".pial \
                    "$SUBJECTS_DIR"/"$subject"/surf/"$x".sphere.reg \
                    "$templateDir"conte69/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemisphere^^}.sphere.32k_fs_LR.surf.gii \
                    "$SUBJECTS_DIR"/"$subject"//surf/"$x".midthickness.surf.gii \
                    "$outDir"/"$subject"_"$x"_midthickness_32k_fs_LR.surf.gii \
                    "$outDir"/"$subject"_"$x"_sphereReg.surf.gii

    wb_command -metric-resample \
        "$tmpProcessingDirectory"/"$subject"_fmri2fs_bbr_"$x".func.gii \
        "$outDir"/"$subject"_"$x"_sphereReg.surf.gii \
        "$templateDir"conte69/resample_fsaverage/fs_LR-deformed_to-fsaverage.${hemisphere^^}.sphere.32k_fs_LR.surf.gii \
        ADAP_BARY_AREA \
        "$tmpProcessingDirectory"/"$subject"_fmri2fs_bbr_"$x"_c69-32k.func.gii \
        -area-surfs \
        "$SUBJECTS_DIR"/"$subject"/surf/"$x".midthickness.surf.gii \
        "$outDir"/"$subject"_"$x"_midthickness_32k_fs_LR.surf.gii
    mri_convert "$tmpProcessingDirectory"/"$subject"_fmri2fs_bbr_"$x"_c69-32k.func.gii "$outDir"/"$subject"_fmri2fs_bbr_"$x"_c69-32k.mgh
done
