load(['/home/richard/Downloads/micasoft/parcellations/fs_LR-conte69/SJH_1015_32k_fs_LR.mat'])

inDir='/Data2/CamCAN/vol2surf_holdout/';
FS7toC69 = dlmread('/Data2/CamCAN/FS7toC69.txt');

holdlist = strsplit(fileread('/Data2/CamCAN/holdout.txt'));

Z = zeros(1012, 1012, (length(holdlist)-1));
parfor s = 1:length(holdlist)-1
    tic
    lh = load_mgh([inDir '/' holdlist{s} '_fmri2fs_bbr_lh_c69-32k.mgh']);
    rh = load_mgh([inDir '/' holdlist{s} '_fmri2fs_bbr_lh_c69-32k.mgh']);

    onC69 = [squeeze(lh)' squeeze(rh)'];

    ts_parc = Bo_surfData2parcelData(onC69, SonG);
    R = corr(ts_parc);
    Z(:,:,s) = tanh(R);

    %dlmwrite(['/lustre/archive/q10027/CamCAN_BackUp/corrmat/' holdlist{s} '_Z.txt'], Z(:,:,s))
    toc
end

avgGroupZ = nanmean(Z,3);
dlmwrite([inDir '/group_average_Z.txt'], avgGroupZ)
normAngle = connectivity2normangle(avgGroupZ, 90);
[groupEmbedding, results] = mica_diffusionEmbedding(normAngle);
dlmwrite([inDir '/group_average_Z.txt'], avgGroupZ)

 f = figure('units','normalized','outerposition',[0 0 1 1]);
    for g  = 1:3
        % group gradient
        ROnSurf = BoSurfStatMakeParcelData(groupEmbedding(:,g), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [(0.3*g)-0.25 0.7 0.15 0.15], [(0.3*g)-0.13 0.7 0.15 0.15], ...
            1, 2, [min(groupEmbedding(:,g)) max(groupEmbedding(:,g))], flipud(cmap_grad));
        hbar = colorbar('horizontal');
        hbar.Position = [(0.3*g)-0.165 0.7 0.1 0.01];
    end

    figure;
    scatter(groupEmbedding(:,1), groupEmbedding(:,2), 20, yeo_parc, 'filled')
    colormap(cmap_yeo)
