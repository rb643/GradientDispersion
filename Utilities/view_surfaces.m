%% view_surfaces.m 
%
%
% written by Casey Paquola, MICA lab, Montreal
% April 2019
%

close all 

% init project 
for init_project = 1 
    GH          = '/data_/mica1/03_projects/casey/';
    fsDir       = '/data_/mica1/01_programs/Freesurfer-6.0/';
    
    % useful scripts
    addpath([GH '/micasoft/matlab/useful'])
    addpath([GH '/micaopen/MPC'])
    addpath([GH '/micaopen/surfstat_addons']) 
    addpath([GH '/micaopen/surfstat_chicago'])
    addpath([GH '/micasoft/matlab/gifti-1.6/'])
    
    % colormaps
    addpath([GH '/micasoft/matlab/colormaps/viridis']) 
    addpath([GH '/micasoft/matlab/colormaps/cbrewer']) 
    addpath(genpath([GH '/micasoft/matlab/diffusionEmbedding']))
    m = load([GH '/micasoft/matlab/colormaps/margulies_gradient.mat']);
    c = load([GH '/micasoft/matlab/colormaps/cbrewer/colorbrewer.mat']);
    margulies   = m.margulies; 
    colorbrewer = c.colorbrewer; 
    inferno     = inferno(); 
    viridis     = viridis(); 
    magma       = magma(); 
    plasma      = plasma();
    fparul      = fake_parula(); 
    
    % surface and parcellation
    FS = SurfStatAvSurf({[fsDir '/subjects/fsaverage/surf/lh.pial'] [fsDir '/subjects/fsaverage/surf/rh.pial']});
    [~, lh_annot_lab, ~] = read_annotation([GH '/micasoft/parcellations/fsaverage7/lh.sjh.annot']);
    [~, rh_annot_lab, ~] = read_annotation([GH '/micasoft/parcellations/fsaverage7/rh.sjh.annot']);
    parc = [lh_annot_lab; rh_annot_lab];
    uparcel = unique(parc);


end

%% LOAD SUBJECTS
inDir='/data_/mica1/03_projects/casey/sandbox1/CamCAN1/';
sub='sub-CC723197';
num_surfs=12;

% load values
clear MP_lh
clear MP_rh
for ii = 1:num_surfs
    MP_lh(ii,:) = load_mgh(char(strcat(inDir, sub, '/surfaces/equivSurfs/14surfs_out/', sub, '.lh.', string(ii), '.mgh')));
    MP_rh(ii,:) = load_mgh(char(strcat(inDir, sub, '/surfaces/equivSurfs/14surfs_out/', sub, '.rh.', string(ii), '.mgh')));
end
MP = [MP_lh MP_rh];

% load surface
S = SurfStatAvSurf({strcat(inDir, sub, '/surfaces/', sub, '/surf/lh.pial'), strcat(inDir, sub, '/surfaces/', sub, '/surf/rh.pial')});

% visualise
MTmoment(1,:) = [MP_lh(10,:) MP_rh(10,:)];
ksdensity(MTmoment)
close all
f = figure;
SurfStatViewData(MTmoment, S)
SurfStatColLim([0 2])
colormap(colorbrewer.div.Spectral{1,9}/255)

% create mpc matrix (and nodal intensity profiles if parcellating)
parc_name='sjh';
[~, lh_parc, ~] = read_annotation(strcat(inDir, sub, '/surfaces/', sub, '/label/lh.', parc_name, '.annot'));
[~, rh_parc, ~] = read_annotation(strcat(inDir, sub, '/surfaces/', sub, '/label/rh.', parc_name, '.annot'));
[MPC, I] = build_mpc(MP, vertcat(lh_parc,rh_parc));

% Project onto average surface
f = figure;
SurfStatViewData(BoSurfStatMakeParcelData(I(1,:), FS, parc), FS)