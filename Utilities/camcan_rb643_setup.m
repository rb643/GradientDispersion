% init project
for init_project = 1
  % annoying issue with multple matlab installations
  restoredefaultpath
  rehash toolboxcache
  savepath

    GH = '/home/richard/Dropbox/Tools';

    %mainDir     = [GH '/sandbox1/CamCAN/'];
    inDir       = '/Data2/CamCAN/vol2surf/';
    gitRepo     = ['/Data2/yesweCanCAM/'];

    % useful scripts
    addpath(genpath(gitRepo));
    addpath([gitRepo '/Utilities/'])
    addpath([GH '/micasoft/matlab/vtk']);
    addpath([GH '/micasoft/matlab/NIfTI_20140122']);
    addpath([GH '/micasoft/matlab/useful'])
    addpath([GH '/micaopen/surfstat_addons'])
    addpath([GH '/micaopen/surfstat_chicago'])
    addpath(genpath([GH '/micasoft/matlab/diffusionEmbedding']))
    addpath([GH '/cifti-matlab'])
    addpath([GH '/gifti-matlab'])
    addpath([GH '/BCT/'])
    addpath([GH '/export_fig/'])
    addpath([GH '/micasoft/sandbox/casey'])
    addpath([GH '/matlab_bgl'])
    addpath('/usr/local/freesurfer/matlab/')

    % colormaps
    addpath([GH '/micasoft/matlab/colormaps/viridis'])
    addpath([GH '/micasoft/matlab/colormaps/cbrewer'])
    m = load([GH '/micasoft/matlab/colormaps/margulies_gradient.mat']);
    c = load([GH '/micasoft/matlab/colormaps/cbrewer/colorbrewer.mat']);
    margulies   = m.margulies;
    colorbrewer = c.colorbrewer;
    inferno     = inferno();
    viridis     = viridis();
    magma       = magma();
    plasma      = plasma();
    fparul      = fake_parula();
    load([GH '/ScientificColourMaps5/bamako/bamako.mat'])
    cmap_grad = bamako;
    load([GH '/ScientificColourMaps5/oslo/oslo.mat'])
    cmap_grad2 = oslo;
    load([GH '/ScientificColourMaps5/lajolla/lajolla.mat'])
    cmap_strength = lajolla;
    cmap_var  = flipud(gray);
    cmap_beta = flipud(interp_colormap(colorbrewer.div.RdBu{1,9}/255, 23));
    load([gitRepo '/Utilities/cmap_yeo.mat'])

    % load surfaces
    load([GH '/micasoft/parcellations/fs_LR-conte69/fsaverage.midthickness_mni_32k_fs_LR.mat'])
    load([GH '/micasoft/parcellations/fs_LR-conte69/SJH_1015_32k_fs_LR.mat'])
    parc = SonG;
    uparcel = unique(parc);

    % wall mask
    load([gitRepo '/Utilities/corticalMask.mat'])
    mask = mask.cortex;

    % atlases
    lh_yeo = gifti([GH '/micasoft/parcellations/fs_LR-conte69/maps/lh.Yeo2011_7Networks_N1000_32k.label.gii']);
    rh_yeo = gifti([GH '/micasoft/parcellations/fs_LR-conte69/maps/rh.Yeo2011_7Networks_N1000_32k.label.gii']);
    yeo = [lh_yeo.cdata; rh_yeo.cdata];
    uparcel = unique(parc);
    yeo_parc = zeros(size(uparcel));
    for node = 1:length(uparcel)
        thisparcel      = uparcel(node);
        yeo_parc(node) = mode(yeo(parc==thisparcel));
    end
end
