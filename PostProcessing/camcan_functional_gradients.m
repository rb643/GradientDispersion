%% camcan_functional_gradients.m
%
%
% written by Casey Paquola, MICA lab, Montreal
% November 2018
%

% init project
for init_project = 1

    GH = '/data_/mica1/03_projects/casey/';

    mainDir     = [GH '/sandbox1/CamCAN/'];
    inDir       = '/data_/mica2/CamCAN/';
    gitRepo     = [GH '/yesweCanCAM/'];

    % useful scripts
    addpath(genpath(mainDir));
    addpath(genpath(gitRepo));
    addpath([GH '/micasoft/matlab/vtk']);
    addpath([GH '/micasoft/matlab/NIfTI_20140122']);
    addpath([GH '/micasoft/matlab/useful'])
    addpath([GH '/micaopen/MPC'])
    addpath([GH '/micaopen/surfstat_addons'])
    addpath([GH '/micaopen/surfstat_chicago'])
    addpath(genpath([GH '/micasoft/matlab/diffusionEmbedding']))
    addpath('/data_/mica1/02_codes/matlab_toolboxes/cifti-matlab')
    addpath('/data_/mica1/02_codes/BCT/')

    % colormaps
    cmap_dir = [GH '/micasoft/matlab/colormaps/'];
    addpath([cmap_dir '/viridis'])
    addpath([cmap_dir '/cbrewer'])
    m = load([cmap_dir '/margulies_gradient.mat']);
    c = load([cmap_dir '/cbrewer/colorbrewer.mat']);
    margulies   = m.margulies;
    colorbrewer = c.colorbrewer;
    inferno     = inferno();
    viridis     = viridis();
    magma       = magma();
    plasma      = plasma();
    fparul      = fake_parula();
    load([cmap_dir '/ScientificColourMaps5/bamako/bamako.mat'])
    cmap_grad = bamako;
    load([cmap_dir 'ScientificColourMaps5/oslo/oslo.mat'])
    cmap_grad2 = oslo;
    load([cmap_dir 'ScientificColourMaps5/lajolla/lajolla.mat'])
    cmap_strength = lajolla;
    cmap_var  = flipud(gray);
    cmap_beta = flipud(interp_colormap(colorbrewer.div.RdBu{1,9}/255, 23));
    load([cmap_dir 'cmap_yeo.mat']);

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


%% Load data
build_gradients = 0;

for import_data = 1

    if build_gradients == 1

        holdOut = dlmread([gitRepo '/Data/group_embedding.txt']);

        subjects = strsplit(fileread('subject_list.txt'));
        for s = 1:length(subjects)


            try
                lh = load_mgh(strcat(inDir, 'sub-', subjects{s}, '_fmri2fs_bbr_lh_c69-32k.mgh'));
                rh = load_mgh(strcat(inDir, 'sub-', subjects{s}, '_fmri2fs_bbr_rh_c69-32k.mgh'));
                ts = [squeeze(lh); squeeze(rh)];

                ts_parc = zeros(256, length(uparcel));
                for ii = 1:length(uparcel)
                    this_parc = uparcel(ii);
                    ts_parc(:,ii) = mean(ts(parc==this_parc,:));
                end

                R = corr(ts_parc);
                Z = tanh(R);
                dlmwrite([inDir '/func_mat/' subjects{s} '_Z.txt'], Z)

                % diffusion embedding
                normAngle = connectivity2normangle(Z, 90);
                [indiEmbed, results] = mica_diffusionEmbedding(normAngle);

                % match with original embedding and store lambdas
                r = corr([indiEmbed(:,1:3) holdOut(:,1:3)]);
                [~, idx] = max(abs(r(4,1:3)));
                indiVar(s,1) = results.lambdas(idx)/sum(results.lambdas);
                [~, idx] = max(abs(r(5,1:3)));
                indiVar(s,2) = results.lambdas(idx)/sum(results.lambdas);
                [~, idx] = max(abs(r(6,1:3)));
                indiVar(s,3) = results.lambdas(idx)/sum(results.lambdas);

                % align
                [U,~,V] = svd(holdOut.' * indiEmbed,0);
                U = U(:, 1:min(size(U,2), size(V,2)));
                V = V(:, 1:min(size(U,2), size(V,2)));
                xfms = V * U';
                alignEmbed(:,:,s) = indiEmbed * xfms;

                % write out
                dlmwrite([inDir '/func_grad/' subjects{s} '_gradient.txt'], alignEmbed(:,:,s))

                sub_successs(s) = 1;
            catch
                sub_successs(s) = 0;
            end

        end


    else

            flist = dir([gitRepo '/Data/gradients/']);
            sublist = cell(1, length(flist)-2);
            G1      = zeros(length(uparcel), 28, length(flist)-2);
            for s  = 1:length(flist)-2
                G1(:,:,s) = dlmread([gitRepo '/Data/gradients/' flist(s+2).name]);
                tmp = strsplit(flist(s+2).name, '_');
                sublist{s} = tmp{1};
            end

    end

    % load demographic variables
    [~, ~, pheno]   = xlsread([gitRepo '/Data/Health.xlsx']);
    mot_list        = dir([gitRepo '/Data/motion_files']);
    age     = zeros(length(sublist), 1);
    sex     = cell(length(sublist), 1);
    motion  = zeros(length(sublist), 1);
    for ii = 1:length(sublist)
        for jj = 1:length(pheno)
            if strcmp(sublist{ii}, pheno{jj,3})
                age(ii)     = pheno{jj,7};
                sex{ii}   = pheno{jj,11};
            end
        end
        motion(ii) = mean(fileread([gitRepo '/Data/motion_files/sub-' sublist{ii} '_task-Rest_bold_metric_REFRMS.1D']));
    end

end

%% Figure 1a
for group_gradients = 1

    f = figure('units','normalized','outerposition',[0 0 1 1]);
    for g  = 1:3

        grad = reshape(G1(:,g,:), [], length(sublist));

        % run pca for group gradient
        group_grad(:,g)  = mean(grad,2);

        % group gradient
        ROnSurf = BoSurfStatMakeParcelData(group_grad(:,g), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [(0.3*g)-0.25 0.7 0.15 0.15], [(0.3*g)-0.13 0.7 0.15 0.15], ...
            1, 2, [min(group_grad(:,g)) max(group_grad(:,g))], flipud(plasma));
        hbar = colorbar('horizontal');
        hbar.Position = [(0.3*g)-0.165 0.7 0.1 0.01];

        % variance in group gradient
        ROnSurf = BoSurfStatMakeParcelData(std(grad, [], 2), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [(0.3*g)-0.25 0.52 0.15 0.15], [(0.3*g)-0.13 0.52 0.15 0.15], ...
            1, 2, [min(std(grad,[],2)) max(std(grad,[],2))], cmap_var);
        hbar = colorbar('horizontal');
        hbar.Position = [(0.3*g)-0.165 0.52 0.1 0.01];

    end
    exportfigbo(f, [gitRepo '/Figures/subplots/Figure_1A.png'], 'png', 10)


end

for node_changes = 1

    f = figure('units','normalized','outerposition',[0 0 1 1]);

    %node wise shifts in gradient values
    nbin = 10;
    for g  = 1:3

        grad = reshape(G1(:,g,:), [], length(sublist));

        % estimate regional effects
        for k = 1:size(G1,1)

            T = table(reshape(G1(k,g,:), length(sublist), []), ...
                age, motion, sex, ...
                'VariableNames', {'grad', 'age', 'motion', 'sex'});

            mdl = fitlm(T, 'grad ~ age + motion + sex');

            for c = 1:4
                mdl_out(k,c) = table2array(mdl.Coefficients(2,c));
            end
            %mdl_out(k,5) = mdl.RMSE;
            %mdl_out(k,6) = mdl.Rsquared.Adjusted;
            %mdl_resid(k,:) = mdl.Residuals.Standardized;

            % write them to a csv files for other analyses and posterity
            T = array2table(mdl_out,...
                'VariableNames', {'lowerCI', 'upperCI', 'T', 'P'});
          writetable(T, char(strcat(gitRepo, '/Data/age_effect_G', num2str(g), '.csv')));

        end

        pcorr = fdr_bh(mdl_out(:,4), 0.0083);

        % projection of beta estimates on the surface
        ROnSurf = BoSurfStatMakeParcelData((mdl_out(:,1).*pcorr), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [(0.3*g)-0.25 0.7 0.15 0.15], [(0.3*g)-0.13 0.7 0.15 0.15], ...
            1, 2, [-(max(abs(mdl_out(:,1)))) max(abs(mdl_out(:,1)))], cmap_beta);
        hbar = colorbar('horizontal');
        hbar.Position = [(0.3*g)-0.165 0.7 0.1 0.01];

        % rank
        [~, idx] = sort(group_grad(:,g));
        Grank = sort_back([1:size(group_grad,1)]', idx);

        %' group-level gradient on the surface
        vGrad = BoSurfStatMakeParcelData(group_grad(:,g), G, parc);
        vs = SurfStatSmooth(vGrad, G, 5);
        extent = [min(vs) max(vs)];

        % estimate age effect within each bin
        for b = 1:nbin

            % select nodes within bin - accounts for equal number of nodes per
            % bin
            idx = ((Grank < round((length(Grank) / nbin) * b))...
                + (Grank > round((length(Grank) / nbin) * (b - 1)))) == 2;

            % visualise bins
            binMask = BoSurfStatMakeParcelData(idx, G, parc);
            Bs = vs.* mask .* binMask;
            Bs(Bs==0) = (min(vs) - 0.01);
            BoSurfStat_calibrate2Views(Bs, G, [(b*0.025)+(g*0.3)-0.26 0.625 0.025 0.025], [(b*0.025)+(g*0.3)-0.26 0.6 0.025 0.025], ...
                1, 2, extent, [1 1 1; cmap_grad])

            % linear model within bin
            binval = reshape(mean(grad(idx,:)), length(sublist), 1);
            T = table(binval, age, motion, sex, ...
                'VariableNames', {'grad', 'age', 'motion', 'sex'});
            mdl = fitlm(T, 'grad ~ age + motion + sex');
            binlm(b) = table2array(mdl.Coefficients(2,1));

        end

        % beta estimates across bins
        a(g) = axes('position', [(g*0.3)-0.22 0.42 (0.0225*nbin) 0.15]);
        plot(1:nbin, binlm, 'k', 'LineWidth', 2); hold on;
        scatter(1:nbin, binlm, 70, binlm, 'filled');
        colormap(a(g), cmap_beta);
        xlim([1 10])

    end
    exportfigbo(f, [gitRepo '/Figures/F1_nodal_changes.png'], 'png', 10)
end


%% Figure 1b - distribution of gradient values
for f1b_statistics = 1
    % calculate metrics
    r1 = zeros(length(sublist), 3);
    h1 = zeros(length(sublist), 3);
    dip = zeros(length(sublist), 3);
    skew = zeros(length(sublist), 3);
    kurt = zeros(length(sublist), 3);
    mdl_hist = zeros(5, 4, 3);
    for g = 1:3

        for s = 1:length(sublist)

            % gradient range
            r1(s,g) = range(G1(:,g,s));

            % create histogram
            [f, xi] = ksdensity(G1(:,g,s));

            % find peaks
            [pks, locs] = findpeaks(f);

            % find height of first and last peak
            h1(s,g)     = double(pks(1));

            % distribution characterisation
            dip(s,g) = HartigansDipTest(f);
            skew(s,g) = skewness(f);
            kurt(s,g) = kurtosis(f);

        end

        T = table(age, sex, motion, r1(:,g), h1(:,g), dip(:,g), skew(:,g), kurt(:,g),...
            'VariableNames', {'age', 'sex', 'motion', 'r1', 'h1', 'dip', 'skew', 'kurt'});
        writetable(T, char(strcat(gitRepo, '/Data/bimodality_values_gradient', num2str(g), '.csv')));

        mdl = fitlm(T, 'r1 ~ age + motion + sex');
        for c = 1:4
            mdl_hist(1,c,g) = table2array(mdl.Coefficients(2,c));
        end
        mdl = fitlm(T, 'h1 ~ age + motion + sex');
        for c = 1:4
            mdl_hist(2,c,g) = table2array(mdl.Coefficients(2,c));
        end
        mdl = fitlm(T, 'dip ~ age + motion + sex');
        for c = 1:4
            mdl_hist(3,c,g) = table2array(mdl.Coefficients(2,c));
        end
            mdl = fitlm(T, 'skew ~ age + motion + sex');
        for c = 1:4
            mdl_hist(4,c,g) = table2array(mdl.Coefficients(2,c));
        end
        mdl = fitlm(T, 'kurt ~ age + motion + sex');
        for c = 1:4
            mdl_hist(5,c,g) = table2array(mdl.Coefficients(2,c));
        end

    end
end

for f1b_joyplotting = 1
    % decile average histograms
    age_decile = (age>=18) + (age>=28) + (age>=38) + (age>=48) + (age>=58) + (age>=68) + (age>=78) + (age>=88);
    for g = 1:3
        for ii = 1:max(age_decile)
            decile_grad(:,ii) = mean(G1(:,g,age_decile==ii),3);
        end
        dlmwrite(char(strcat(gitRepo, '/Data/deciles_gradient', num2str(g), '.txt')), decile_grad)
        end
end






%% Figure 2 - age-related changes in yeo networks across 3d space

for f2_analysis = 1

    % find centre of gravity of each yeo network in 3D space
    yeo_cog     = zeros(length(sublist), 3, 7);
    yeo_netSSD  = zeros(length(sublist), 7);
    mdl_yeo     = zeros(2, 4, 7);
    for n = 1:7
        netNodes = find(yeo_parc==n);
        group_yeo_cog(n,:) = median(group_grad(netNodes,:));
        for s = 1:length(sublist)
            % centre of gravity
            yeo_cog(s,:,n) = median(G1(netNodes, 1:3, s));
            % distance from group average
            yeo_cog_dist(s,n) = pdist2(group_yeo_cog(n,:), yeo_cog(s,:,n));
            % sum squared distance of nodes to network cog - within network
            % integration
            clear dist_to_centre
            for ii = 1:length(netNodes)
                dist_to_centre(ii) = pdist2(G1(netNodes(ii),1:3,s), yeo_cog(s,:,n));
            end
            yeo_netSSD(s,n) = sumsqr(dist_to_centre);
        end


        T = table(age, sex, motion, yeo_cog(:,1,n), yeo_cog(:,2,n), yeo_cog(:,3,n), ...
            yeo_cog_dist(:,n), 1./yeo_netSSD(:,n), 'VariableNames', ...
            {'age', 'sex', 'motion', 'g1', 'g2', 'g3', 'dist', 'disp'});
        mdl = fitlm(T, 'g1 ~ age + motion + sex');
        for c = 1:4
            mdl_yeo(1,c,n) = table2array(mdl.Coefficients(2,c));
        end
        mdl = fitlm(T, 'g2 ~ age + motion + sex');
        for c = 1:4
            mdl_yeo(2,c,n) = table2array(mdl.Coefficients(2,c));
        end
        mdl = fitlm(T, 'g3 ~ age + motion + sex');
        for c = 1:4
            mdl_yeo(3,c,n) = table2array(mdl.Coefficients(2,c));
        end
        mdl = fitlm(T, 'disp ~ age + motion + sex');
        for c = 1:4
            mdl_yeo(4,c,n) = table2array(mdl.Coefficients(2,c));
        end
    end

    % between network segregation
    net_dist = zeros(7,7,length(sublist));
    net_overlap = zeros(7,7,length(sublist));
    for n1 = 1:7
        for n2 = 1:7
            % distance between network cog
            net_dist(n1,n2,:) = sqrt(sum((yeo_cog(:,:,n1)' - yeo_cog(:,:,n2)').^2));
            % overlap of networks
            n1Nodes = find(yeo_parc==n1);
            n2Nodes = find(yeo_parc==n2);
            for s = 1:length(sublist)
                these_nodes = G1([n1Nodes n2Nodes],1:3,s);
                idx = kmeans(these_nodes,2);
                yeo_idx = [ones(length(n1Nodes),1); ones(length(n2Nodes),1)*2];
                net_overlap(n1,n2,s) = pdist([idx'; yeo_idx'], 'hamming');
            end
        end
    end
    for s = 1:length(sublist)
        net_dist_long(s,:) = squareform(net_dist(:,:,s));
        tmp = net_overlap(:,:,s);
        tmp = tmp .* ~eye(size(tmp));
        net_overlap_long(s,:) = squareform(tmp);
    end
    for ii = 1:21
        T = table(age, sex, motion, net_dist_long(:,ii), net_overlap_long(:,ii), ...
            'VariableNames', {'age', 'sex', 'motion', 'dist', 'over'});
        mdl = fitlm(T, 'dist ~ age + sex + motion');
        dist_tstat(ii) = table2array(mdl.Coefficients(2,3));

        mdl = fitlm(T, 'over ~ age + sex + motion');
        over_tstat(ii) = table2array(mdl.Coefficients(2,3));
    end

end

for f2_statistics = 1

  % null model for change in SSD
  load([gitRepo '/Utilities/perm_sphere_10000.mat']);
  n_perm =  1000;
  perm_yeo_disp = zeros(n_perm, 7);
  perm_dist_tstat = zeros(n_perm, 21);
  perm_over_tstat = zeros(n_perm, 21);
  clear perm_dist_tstat
  clear perm_over_tstat
  for create_null = 1
    parfor_progress(n_perm)
    parfor p = 1:n_perm
      yeo_rot = yeo_parc(perm_all(:,p));
      perm_yeo_cog    = zeros(length(sublist), 3, 7);
      perm_yeo_netSSD = zeros(length(sublist), 7);
      for n = 1:7
        netNodes = find(yeo_rot==n);

        % sum squared distance of each network
        for s = 1:length(sublist)
          % centre of gravity
          perm_yeo_cog(s,:,n) = median(G1(netNodes, 1:3, s));
          % sum squared distance of nodes to network cog - within network
          % integration
          dist_to_centre = zeros(1, length(netNodes));
          for ii = 1:length(netNodes)
            dist_to_centre(ii) = pdist2(G1(netNodes(ii),1:3,s), perm_yeo_cog(s,:,n));
          end
          perm_yeo_netSSD(s,n) = sumsqr(dist_to_centre);
        end
        T = table(age, sex, motion, 1./perm_yeo_netSSD(:,n), ...
        'VariableNames', {'age', 'sex', 'motion', 'disp'});
        mdl = fitlm(T, 'disp ~ age + motion + sex');
        perm_yeo_disp(p,n) = table2array(mdl.Coefficients(2,3));
      end

      % between network segregation
      perm_net_dist = zeros(7,7,length(sublist));
      perm_net_over = zeros(7,7,length(sublist));
      for n1 = 1:7
        for n2 = 1:7
          % distance between network cog
          perm_net_dist(n1,n2,:) = sqrt(sum((perm_yeo_cog(:,:,n1)' - perm_yeo_cog(:,:,n2)').^2));
          % overlap of networks
          n1Nodes = find(yeo_rot==n1);
          n2Nodes = find(yeo_rot==n2);
          for s = 1:length(sublist)
            these_nodes = G1([n1Nodes n2Nodes],1:3,s);
            idx = kmeans(these_nodes,2);
            yeo_idx = [ones(length(n1Nodes),1); ones(length(n2Nodes),1)*2];
            perm_net_over(n1,n2,s) = pdist([idx'; yeo_idx'], 'hamming');
          end
        end
      end
      perm_net_dist_long      = zeros(length(sublist), 21);
      perm_net_overlap_long   = zeros(length(sublist), 21);
      for s = 1:length(sublist)
        perm_net_dist_long(s,:) = squareform(perm_net_dist(:,:,s));
        tmp = perm_net_over(:,:,s);
        tmp = tmp .* ~eye(size(tmp));
        perm_net_overlap_long(s,:) = squareform(tmp);
      end
      for ii = 1:21
        T = table(age, sex, motion, perm_net_dist_long(:,ii), perm_net_overlap_long(:,ii), ...
        'VariableNames', {'age', 'sex', 'motion', 'dist', 'over'});
        mdl = fitlm(T, 'dist ~ age + sex + motion');
        perm_dist_tstat{p}(ii) = table2array(mdl.Coefficients(2,3));
        mdl = fitlm(T, 'over ~ age + sex + motion');
        perm_over_tstat{p}(ii) = table2array(mdl.Coefficients(2,3));
      end
     parfor_progress;
    end %end parfor
  end

  % save out results for integration (within network SSD)
  yeo_SSD_age_p = sum(perm_yeo_disp>(squeeze(mdl_yeo(4,3,:))'))/n_perm;
  dlmwrite([gitRepo '/Data/yeo_SSD_age.csv'], [ squeeze(mdl_yeo(4,3,:))'; perm_yeo_disp])
  dlmwrite([gitRepo '/Data/yeo_SSD_age_p.csv'], yeo_SSD_age_p)

  % save out results for segregation (distance between and overlap of
  % networks)
  yeo_dist_age_p = sum(cell2mat(perm_dist_tstat')>dist_tstat)/n_perm;
  dlmwrite([gitRepo '/Data/yeo_dist_age.csv'], [dist_tstat; cell2mat(perm_dist_tstat')])
  dlmwrite([gitRepo '/Data/yeo_dist_age_p.csv'], yeo_dist_age_p)
  yeo_over_age_p = sum(cell2mat(perm_over_tstat')>over_tstat)/n_perm;
  dlmwrite([gitRepo '/Data/yeo_over_age.csv'], [over_tstat; cell2mat(perm_over_tstat')])
  dlmwrite([gitRepo '/Data/yeo_over_age_p.csv'], yeo_over_age_p)

  for create_null = 0
      perm_dist_tstat = csvread([gitRepo '/Data/yeo_dist_age.csv']);
      dist_tstat = squeeze(perm_dist_tstat(1,:));
      perm_dist_tstat = squeeze(perm_dist_tstat(2:end,:));
      yeo_dist_age_p = squeeze(csvread([gitRepo '/Data/yeo_dist_age_p.csv']));
      perm_over_tstat = csvread([gitRepo '/Data/yeo_over_age.csv']);
      over_tstat = squeeze(perm_over_tstat(1,:));
      perm_over_tstat = squeeze(perm_over_tstat(2:end,:));
      yeo_over_age_p = squeeze(csvread([gitRepo '/Data/yeo_over_age_p.csv']));
  end

end

for f2_visualisation = 1

    f = figure('units','normalized','outerposition',[0 0 1 1]);

    % 3d scatter plot, coloured by yeo
    a(1) = axes('position', [0.05 0.3 0.4 0.4]);
    scatter3(group_grad(:,1), group_grad(:,2), group_grad(:,3), 50, yeo_parc, 'filled');
    colormap(a(1), cmap_yeo)

    % plot change in 2d space
    combos = combnk(1:3,2);
    for c = 1:length(combos)
        a(c+1) = axes('position', [0.5 0.72-(c*0.14) 0.08 0.1]);
        x = group_grad(:,combos(c,1));
        y = group_grad(:,combos(c,2));
        scatter(x, y, 10, yeo_parc, 'filled', 'MarkerFaceAlpha', 0.2); hold on;
        xlim([min(x) max(x)])
        ylim([min(y) max(y)])
        colormap(a(c+1), cmap_yeo)
        scatter(group_yeo_cog(:,combos(c,1)), group_yeo_cog(:,combos(c,2)), 80, cmap_yeo(2:end,:), 'filled');
        for ii = 1:7
            quiver(group_yeo_cog(ii,combos(c,1)), group_yeo_cog(ii,combos(c,2)), ...
                mdl_yeo(combos(c,1),3,ii)*0.01, mdl_yeo(combos(c,2),3,ii)*0.01, ...
                'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 1, 'AutoScale', 'off')
        end
        hold off;
    end

    % % plot change in 3d space
    % scatter3(group_grad(:,1), group_grad(:,2), group_grad(:,3), 10, yeo_parc, 'filled', 'MarkerFaceAlpha', 0.2); hold on;
    % colormap(cmap_yeo)
    % scatter3(group_yeo_cog(:,1), group_yeo_cog(:,2), group_yeo_cog(:,3), 100, cmap_yeo(2:end,:), 'filled');
    % for ii = 1:7
    %     quiver3(group_yeo_cog(ii,1), group_yeo_cog(ii,2), group_yeo_cog(ii,3), ...
    %         mdl_yeo(1,3,ii)*0.01,mdl_yeo(2,3,ii)*0.01,mdl_yeo(3,3,ii)*0.01, ...
    %         'Color', cmap_yeo(ii+1,:), 'LineWidth', 3)
    % end

    % within network integration
    a(5) = axes('position', [0.62 0.51 0.1 0.15]);
    b = bar(squeeze(mdl_yeo(4,3,:)));
    colormap(a(5), cmap_yeo(2:end,:))

    % segregation
    a(6) = axes('position', [0.62 0.31 0.1 0.15]);
    imagesc(squareform(dist_tstat))
    caxis([-4 4])
    colormap(a(6), cmap_beta)
    a(7) = axes('position', [0.62 0.29 0.1 0.01])
    imagesc([1:7])
    colormap(a(7), cmap_yeo(2:end,))

    exportfigbo(f, [gitRepo '/Figures/F2_yeo_moves.png'], 'png', 10)
    end


%% Figure 3 - subcortical and structural influences
gen_sub_cort_Z = 0;

for load_data = 1

    % l_thal, l_caud, l_put, l_pall, l_hipp, l_amy, l_acc
    % r_thal, r_caud, r_put, r_pall, r_hipp, r_amy, r_acc
    first_id = [10 11 12 13 17 18 26 49 50 51 52 53 54 58];
    sub_cort_Z = zeros(1026, 14, length(sublist));

    parfor_progress(length(sublist))
    parfor s = 1:length(sublist)
        if gen_sub_cort_Z == 1

          subject = sublist{s};
          [Z R] = camcan_sub_cort_Z(subject, GH, inDir, gitRepo);
          sub_cort = Z(:,(end-13):end);
          sub_cort_r = R(:,(end-13):end);
          %dlmwrite([gitRepo '/Data/subcortex_Z/sub-'  subject '_Z.csv'],sub_cort);
          %dlmwrite([gitRepo '/Data/subcortex_R/sub-'  subject '_R.csv'],sub_cort_r);
          sub_cort_Z(:,:,s) = sub_cort;
          sub_cort_R(:,:,s) = sub_cort_r;
        else
            sub_cort_Z(:,:,s) = csvread([gitRepo '/Data/subcortex_Z/sub-'  sublist{s} '_Z.csv']);
        end
    parfor_progress;
    end

end

% spatial overlap of subcortical-cortical connectivity maps with
% ipsilateral gradients
lh_nodes = 1:504;
rh_nodes = 505:1012;
subcortex_num = length(first_id)/2;
for ii = 1:subcortex_num
    for s = 1:length(sublist)
        %sub_corr(:,ii,s) = [sub_cort_Z(lh_nodes,ii,s); sub_cort_Z(rh_nodes,ii+subcortex_num,s)];
        sub_corr(:,ii,s) = [sub_cort_R(lh_nodes,ii,s); sub_cort_R(rh_nodes,ii+subcortex_num,s)];
        for g = 1:3
            sub_grad(s,g,ii) = corr(squeeze(sub_corr(:,ii,s)), squeeze(G1(:,g,s)));
        end
    end
    T = table(age, sex, motion, sub_grad(:,1,ii), sub_grad(:,2,ii), sub_grad(:,3,ii), ...
        'VariableNames', {'age', 'sex', 'motion', 'corr_g1', 'corr_g2', 'corr_g3'});
    mdl = fitlm(T, 'corr_g1 ~ age + motion + sex');
    for c = 1:4
        mdl_sub_corr(1,c,ii) = table2array(mdl.Coefficients(2,c));
    end
    mdl = fitlm(T, 'corr_g2 ~ age + motion + sex');
    for c = 1:4
        mdl_sub_corr(2,c,ii) = table2array(mdl.Coefficients(2,c));
    end
    mdl = fitlm(T, 'corr_g3 ~ age + motion + sex');
    for c = 1:4
        mdl_sub_corr(3,c,ii) = table2array(mdl.Coefficients(2,c));
    end
end

  % edge-wise linear Model
  tmat = nan(size(sub_corr,1),size(sub_corr,2));
  pmat = nan(size(sub_corr,1),size(sub_corr,2));

  for e = 1:size(sub_corr,1)
      for ee = 1:size(sub_corr,2)
          T = table(age, sex, motion, squeeze(sub_corr(e,ee,:)), ...
          'VariableNames', {'age', 'sex', 'motion', 'conn'});
          mdl = fitlm(T, 'conn ~ age + motion + sex');

          tmat(e,ee,:) = table2array(mdl.Coefficients(2,3));
          pmat(e,ee,:) = table2array(mdl.Coefficients(2,4));

      end
  end

  % generate some csv files for plotting
  % group_sub_corr = mean(sub_corr,3);
  net_label = repmat(yeo_parc, 1, 7);
  sub_label = repmat(1:7, 1012, 1);

  out(:,1) = reshape(pmat, size(pmat,1)*size(pmat,2), 1);
  out(:,2) = reshape(net_label, size(pmat,1)*size(pmat,2), 1);
  out(:,3) = reshape(sub_label, size(pmat,1)*size(pmat,2), 1);
  dlmwrite([gitRepo 'Data/subcort_cort_conn_pmat.csv'],out);

% age-related changes in variance of subcortical-yeo connectivity with integration of
% yeo networks
mdl_sub_yeo_med_int = zeros(4,4,7,7);
for n = 1:7
    netNodes = (yeo_parc==n);
    % connectivity to yeo
    sub_yeo_corr(:,:,n) = squeeze(median(sub_corr(netNodes,:,:)))'; % subject rows, subcortical as columns
    for ii = 1:subcortex_num
        T = table(age, sex, motion, sub_yeo_corr(:,ii,n), 1./yeo_netSSD(:,n), ...
            'VariableNames', {'age', 'sex', 'motion', 'sub_yeo_corr', 'disp'});
        % a - age to mediator
        mdl = fitlm(T, 'sub_yeo_corr ~ age + motion + sex');
        for c = 1:4
            mdl_sub_yeo_med_int(1,c,n,ii) = table2array(mdl.Coefficients(2,c));
        end

        % c - age to yeo network integration
        mdl = fitlm(T, 'disp ~ age + motion + sex');
        for c = 1:4
            mdl_sub_yeo_med_int(3,c,n,ii) = table2array(mdl.Coefficients(2,c));
        end

        % b - mediator to yeo network integration
        % c' - age to yeo network integration, controlling for mediator
        mdl = fitlm(T, 'disp ~ age + motion + sex + sub_yeo_corr');
        for c = 1:4
            mdl_sub_yeo_med_int(2,c,n,ii) = table2array(mdl.Coefficients(5,c));
            mdl_sub_yeo_med_int(4,c,n,ii) = table2array(mdl.Coefficients(2,c));
        end

    end
end
% find key effects
clear med_type
for n = 1:7
    for ii = 1:7
          if sum(mdl_sub_yeo_med_int(1:3,4,n,ii) < 0.01)==3
              if abs(mdl_sub_yeo_med_int(3,3,n,ii)) > abs(mdl_sub_yeo_med_int(4,3,n,ii))
                med_type(n,ii) = 2;
              else
                med_type(n,ii) = 1;
              end
          else
              med_type(n,ii) = 0;
          end
    end
end

for f3_visualisation = 1

    f = figure('units','normalized','outerposition',[0 0 1 1]);

    ii = 1; % thalamus
    n = 7;  % dan = 3, van = 4, fp = 6, dmn = 7

    Age = term(age);
    sex1 = double(char(sex));
    Sex = term(sex1);
    Motion = term(motion);
    subcort = sub_yeo_corr(:,ii,n);
    Subcort = term(subcort);
    yeo_ssd = 1./yeo_netSSD(:,n);

    % path a
    a(1) = axes('position', [0.1 0.3 0.2 0.2]);
    inpic=~((isoutlier(subcort)) | (isoutlier(age)));
    scatter(age(inpic), subcort(inpic), 10, 'k', 'filled'); hold on;
    plot_linear_regression(age(inpic), subcort(inpic));
    hold off;

    % path b
    a(2) = axes('position', [0.7 0.3 0.2 0.2]);
    Model = 1 + Age + Sex + Motion;
    slm = SurfStatLinMod(yeo_ssd, Model);
    slm = SurfStatT( slm, Age);
    yeo_corr = (yeo_ssd - slm.X*slm.coef)';
    inpic=~((isoutlier(subcort)) | (isoutlier(yeo_corr')));
    scatter(subcort(inpic), yeo_corr(inpic), 10, 'k', 'filled');  hold on;
    plot_linear_regression(subcort(inpic), yeo_corr(inpic)');
    labeltext = ['T = ' num2str(slm.t) ' -- Effect =' num2str(slm.ef)];
    text(min(subcort), min(yeo_corr), labeltext)
    hold off

    % path c
    a(3) = axes('position', [0.4 0.5 0.2 0.2]);
    inpic=~((isoutlier(age)) | (isoutlier(yeo_ssd)));
    scatter(age(inpic), yeo_ssd(inpic), 10, 'k', 'filled'); hold on
    plot_linear_regression(age(inpic), yeo_ssd(inpic));
    hold off

    % path c'
    a(4) = axes('position', [0.4 0.8 0.2 0.2]);
    Model = 1 + Subcort + Sex + Motion;
    slm = SurfStatLinMod(yeo_ssd, Model);
    slm = SurfStatT( slm, Subcort);
    yeo_corr = (yeo_ssd - slm.X*slm.coef)';
    inpic=~((isoutlier(age)) | (isoutlier(yeo_corr')));
    scatter(age(inpic), yeo_corr(inpic), 10, 'k', 'filled'); hold on;
    plot_linear_regression(age(inpic), yeo_corr(inpic)');
    labeltext = ['T = ' num2str(slm.t) ' -- Effect =' num2str(slm.ef)];
    text(min(age), min(yeo_corr), labeltext)
    hold off

    % subcortical to network map
    this_subc = mean(sub_corr(:,ii,:),3).*(yeo_parc'==n);
    tmp = this_subc;
    tmp(this_subc==0) = [];
    SOnSurf = BoSurfStatMakeParcelData(this_subc, G, parc);
    BoSurfStat_calibrate2Views(SOnSurf, G, [0.4 0.7 0.12 0.12], ...
        [0.5 0.7 0.12 0.12], 1, 2, [min(tmp) max(tmp)], ...
        [0.5 0.5 0.5; cmap_strength])

    % show distribution of network nodes in 3d space
    a(5) = axes('position', [0.85 0.25 0.1 0.1]);
    scatter3(group_grad(:,1), group_grad(:,2), group_grad(:,3), 10, yeo_parc==n, 'filled', 'MarkerFaceAlpha', 0.2); hold on
    colormap(a(5), [0.5 0.5 0.5; cmap_yeo(n+1,:));
    export_fig 'Figures/Figure_3B1.pdf' -transparent
end


%% edge-wise node changes
for edge_changes == 1

      % load connectivity matrices
          flist = dir([gitRepo '/Data/connmat_gsr/']);
          sublist = cell(1, length(flist)-2);
          Z   = zeros(length(uparcel), length(uparcel), length(flist)-2);

          parfor_progress(length(flist)-2)

          parfor s  = 1:length(flist)-2
              Z(:,:,s) = dlmread([gitRepo '/Data/connmat_gsr' '/' flist(s+2).name]);
              tmp = strsplit(flist(s+2).name, '_');
              sublist{s} = tmp{1};
              parfor_progress;
          end
          parfor_progress(0);

      % load demographic variables
      [~, ~, pheno]   = xlsread([gitRepo '/Data/Health.xlsx']);
      mot_list        = dir([gitRepo '/Data/motion_files']);
      age     = zeros(length(sublist), 1);
      sex     = cell(length(sublist), 1);
      motion  = zeros(length(sublist), 1);
      for ii = 1:length(sublist)
          for jj = 1:length(pheno)
              if strcmp(sublist{ii}, pheno{jj,3})
                  age(ii)     = pheno{jj,7};
                  sex{ii}   = pheno{jj,11};
              end
          end
          motion(ii) = mean(fileread([gitRepo '/Data/motion_files/sub-' sublist{ii} '_task-Rest_bold_metric_REFRMS.1D']));
      end


      parfor n = 1:length(sublist)
          mask = tril(true(size(Z(:,:,n))),-1);
          tmp = squeeze(Z(:,:,n));
          trilledmat(n,:) = tmp(mask);
      end

      % run models
      parfor_progress(length(trilledmat))
      parfor i = 1:length(trilledmat)
          %fprintf('\b|\n');
          T = table(trilledmat(:,i), ...
              age, motion, sex, ...
              'VariableNames', {'edge', 'age', 'motion', 'sex'});
          mdl = fitlm(T, 'edge ~ age + motion + sex');
          mdl_out_p(i,:) = table2array(mdl.Coefficients(2,4));
          mdl_out_t(i,:) = table2array(mdl.Coefficients(2,3));
          parfor_progress;
      end
      parfor_progress(0);

      [h, crit_p, adj_ci_cvrg, fdr] = fdr_bh(mdl_out_p,0.001,'pdep','no');

      pmat = nan(length(uparcel),length(uparcel));
      pmat(mask) = mdl_out_p;

      tmat = nan(length(uparcel),length(uparcel));
      tmat(mask) = mdl_out_t;

      pfdrmat = nan(length(uparcel),length(uparcel));
      pfdrmat(mask) = fdr;

      indices = pfdrmat > 0.01;
      pfdrmat(indices) = 1;

      tmat = tril(tmat)+tril(tmat)';
      pfdrmat = tril(pfdrmat)+tril(pfdrmat)';
      logpfdrmat = log(pfdrmat);

      figure;
      [~, idx] = sort (yeo_parc);
      A(1) = axes('position', [0.2 0.2 0.5 0.5]);
      imagesc(pfdrmat(idx,idx)); cb = colorbar; %caxis([min(min(tmat)) abs(min(min(tmat)))]);
      set(gca,'YTickLabel',[]);
      a=get(cb); %gets properties of colorbar
      a.Position; %gets the positon and size of the color bar
      set(cb,'Position',[0.75 0.2 0.02 0.5]);% To change size

      A(2) = axes('position', [0.18 0.2 0.02 0.5]);
      yeo_sort = yeo_parc(idx);
      imagesc(yeo_sort');
      set(gca,'YTickLabel',[]);
      set(gca,'XTickLabel',[]);

      A(3) = axes('position', [0.2 0.7 0.5 0.02]);
      yeo_sort = yeo_parc(idx);
      imagesc(yeo_sort);
      set(gca,'YTickLabel',[]);
      set(gca,'XTickLabel',[]);

      colormap(A(1), (gray));
      colormap(A(2), cmap_yeo) ;
      colormap(A(3), cmap_yeo) ;

      hFig = figure(2);
      set(hFig, 'Position', [0 0 1000 1000])

      exportfigbo(hFig, [gitRepo '/Figures/subplots/Edge_pfdr_meanFC.png'], 'png', 10)

      csvwrite([gitRepo '/Data/edge_t_meanFC.csv'],tmat);
      csvwrite([gitRepo '/Data/edge_p_meanFC.csv'],pmat);
      csvwrite([gitRepo '/Data/edge_p_fdr_meanFC.csv'],pfdrmat);

      % compute within yeo connc
      % null model for change in within and between connectivity
      load([gitRepo '/Utilities/perm_sphere_10000.mat']);
      n_perm =  1000;
      perm_between_tstat = zeros(n_perm, 49);
      perm_within_tstat = zeros(n_perm, 7);

      parfor_progress((n_perm))
      for p = 1:n_perm
          yeo_rot = yeo_parc(perm_all(:,p));
          yeo_within  = zeros(size(Z,3),7);
          yeo_between  = zeros(size(Z,3),7,7);
          for y = 1:7
              for m = 1:size(Z,3)
                  masky = (yeo_rot == y);
                  tmp = Z(masky,masky,m);
                  yeo_within(m,y) = mean(tmp(:));
                  for y2 = 1:7
                      masky2 = (yeo_rot == y2);
                      temp = Z(masky,masky2,m);
                      yeo_between(m,y,y2) = mean(temp(:));
                  end

              end
              T = table(age, sex, motion, squeeze(yeo_within(:,y)), ...
                  'VariableNames', {'age', 'sex', 'motion', 'disp'});
              mdl = fitlm(T, 'disp ~ age + motion + sex');
              perm_within_tstat(p,y) = table2array(mdl.Coefficients(2,3));
          end

          yeo_between = reshape(yeo_between,514,49);

          for b = 1:49
              T2 = table(age, sex, motion, squeeze(yeo_between(:,b)), ...
                  'VariableNames', {'age', 'sex', 'motion', 'disp'});
              mdl = fitlm(T2, 'disp ~ age + motion + sex');
              perm_between_tstat(p,b) = table2array(mdl.Coefficients(2,3));
          end

          display(p);
          %parfor_progress;
      end %end parfor

      csvwrite([gitRepo '/Data/yeo_perm_t_within.csv'],perm_within_tstat);
      csvwrite([gitRepo '/Data/yeo_perm_t_between.csv'],perm_between_tstat);


      yeo_within  = zeros(size(Z,3),7);
      yeo_between  = zeros(size(Z,3),7,7);
      %parfor_progress(size(Z,3))
      for y = 1:7
          for m = 1:size(Z,3)
              masky = (yeo_parc == y);
              tmp = Z(masky,masky,m);
              yeo_within(m,y) = mean(tmp(:));
              for y2 = 1:7
                 masky2 = (yeo_parc == y2);
                 temp = Z(masky,masky2,m);
                 yeo_between(m,y,y2) = mean(temp(:));
              end
        end
        T = table(age, sex, motion, squeeze(yeo_within(:,y)), ...
        'VariableNames', {'age', 'sex', 'motion', 'disp'});
        mdl = fitlm(T, 'disp ~ age + motion + sex');
        within_tstat(y) = table2array(mdl.Coefficients(2,3));

      end
      yeo_between = reshape(yeo_between,515,49);

     for b = 1:49
                T = table(age, sex, motion, squeeze(yeo_between(:,b)), ...
            'VariableNames', {'age', 'sex', 'motion', 'disp'});
        mdl = fitlm(T, 'disp ~ age + motion + sex');
        between_tstat(b) = table2array(mdl.Coefficients(2,3));
        between_pstat(b) = numel(find(between_tstat(:,b)>perm_between_tstat(:,b)))./1000;

    end

      csvwrite([gitRepo '/Data/yeo_within.csv'],yeo_within);
      csvwrite([gitRepo '/Data/yeo_between.csv'],yeo_between);
      csvwrite([gitRepo '/Data/yeo_within_t.csv'],within_tstat);
      csvwrite([gitRepo '/Data/yeo_between_t.csv'],between_tstat);
       csvwrite([gitRepo '/Data/yeo_between_p.csv'],between_pstat);


      % compute graph metrics
      parfor_progress(size(Z,3))
      parfor m = 1:size(Z,3)
          matrix = squeeze(Z(:,:,m));
          matrix=(matrix+matrix')/2;      % force symmetrize matrix
          matrix(matrix<=0) = 0;
           % create random matrices
          rand_matrix = randomize_matrix(matrix);
%           try
%             matrix = mst_threshold(matrix,0.1,'false');
%           catch
%             display('error with the MST thresholding')
%             matrix = nan(1012,1012);
%           end
          Cp = clustering_coef_matrix(matrix, 'O');
          rand_Cp = clustering_coef_matrix(rand_matrix, 'O');

          Gmatrix = sparse(matrix);
          rand_Gmatrix = sparse(rand_matrix);
          dist_matrix = graphallshortestpaths(Gmatrix);
          Lp = mean(dist_matrix);
          rand_dist_matrix = graphallshortestpaths(rand_Gmatrix);
          rand_Lp = mean(rand_dist_matrix);

          % normalised topology
          Cnorm(m,:) = Cp ./ rand_Cp;
          Lnorm(m,:) = Lp ./ rand_Lp;
          parfor_progress;
      end
      parfor_progress(0);

      % check for failures
      outlierL = find(isinf(mean(Lnorm')));
      outlierC = find(isinf(mean(Cnorm')));
      csvwrite([gitRepo '/Data/failed_graph_length.csv'],outlierL);
      csvwrite([gitRepo '/Data/failed_graph_clustering.csv'],outlierC);
      % and remove those
      Lnorm(outlierL,:) = [];
      Cnorm(outlierC,:) = [];

      if subplots == 1
        f = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
        ROnSurf = BoSurfStatMakeParcelData(nanmean(Cnorm), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.1 0.33 0.45 0.45], [0.5 0.33 0.45 0.45], ...
            1, 2, [min(nanmean(Cnorm)) max(nanmean(Cnorm))], flipud(plasma));
        hbar = colorbar('horizontal');
        hbar.Position = [0.3 0.1 0.4 0.1];

        exportfigbo(f, [gitRepo '/Figures/subplots/ClusterMean.png'], 'png', 10)
        csvwrite([gitRepo '/Data/clustering_coeff.csv'],Cnorm);
        close(f)

        exportfigbo(f, [gitRepo '/Figures/subplots/ClusterMean.png'], 'png', 10)
        csvwrite([gitRepo '/Data/clustering_coeff.csv'],Cnorm);
        close(f)


        f = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
        ROnSurf = BoSurfStatMakeParcelData(nanmean(Lnorm), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.1 0.33 0.45 0.45], [0.5 0.33 0.45 0.45], ...
            1, 2, [min(nanmean(Lnorm)) max(nanmean(Lnorm))], flipud(plasma));
        hbar = colorbar('horizontal');
        hbar.Position = [0.3 0.1 0.4 0.1];

        exportfigbo(f, [gitRepo '/Figures/subplots/PathLengthMean.png'], 'png', 10)
        csvwrite([gitRepo '/Data/path_length.csv'],Lnorm);
        close(f)
      end

      %do some age related analyses on clustering
      ageC = age; ageC(outlierC,:) = [];
      sexC = sex; sexC(outlierC,:) = [];
      motionC = motion; motionC(outlierC,:) = [];
      parfor_progress(length(Cnorm))
      parfor i = 1:length(Cnorm)
          T = table(Cnorm(:,i), ...
              ageC, motionC, sexC, ...
              'VariableNames', {'edge', 'age', 'motion', 'sex'});
          mdl = fitlm(T, 'edge ~ age + motion + sex');
          mdl_c_p(i,:) = table2array(mdl.Coefficients(2,4));
          mdl_c_t(i,:) = table2array(mdl.Coefficients(2,3));
          parfor_progress;
      end
      parfor_progress(0);
      csvwrite([gitRepo '/Data/mdl_c_t.csv'],mdl_c_t);
      csvwrite([gitRepo '/Data/mdl_c_p.csv'],mdl_c_p);


      [h, crit_p, adj_ci_cvrg, mdl_c_pfdr] = fdr_bh(mdl_c_p,0.001,'pdep','no');
      mdl_c_t(mdl_c_pfdr > 0.01) = 0;
      mdl_c_pfdr(mdl_c_pfdr > 0.01) = 1;

      if subplots == 1

      % plot
      f = figure('units','normalized','outerposition',[0 0 0.55 0.5]);

        ROnSurf = BoSurfStatMakeParcelData(mdl_c_t, G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.0 0.33 0.30 0.30], [0.25 0.33 0.30 0.30], ...
            1, 2, [-max(mdl_c_t) max(mdl_c_t)], cmap_beta);
        hbar = colorbar('horizontal');
        hbar.Position = [0.12 0.1 0.3 0.1];

        ROnSurf = BoSurfStatMakeParcelData((mdl_c_pfdr), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.5 0.33 0.30 0.30], [0.75 0.33 0.30 0.30], ...
            1, 2, [0 1], plasma);
        hbar = colorbar('horizontal');
        hbar.Position = [0.62 0.1 0.3 0.1];
        exportfigbo(f, [gitRepo '/Figures/subplots/AgeChange_Clustering.png'], 'png', 10)
        close(f)
      end


      %% age models for path length
      ageL = age; ageL(outlierL,:) = [];
      sexL = sex; sexL(outlierL,:) = [];
      motionL = motion; motionL(outlierL,:) = [];
      parfor_progress(length(Lnorm))
      parfor i = 1:length(Lnorm)
          T = table(Lnorm(:,i), ...
              ageL, motionL, sexL, ...
              'VariableNames', {'edge', 'age', 'motion', 'sex'});
          mdl = fitlm(T, 'edge ~ age + motion + sex');
          mdl_l_p(i,:) = table2array(mdl.Coefficients(2,4));
          mdl_l_t(i,:) = table2array(mdl.Coefficients(2,3));
          parfor_progress;
      end
      parfor_progress(0);

      csvwrite([gitRepo '/Data/mdl_l_t.csv'],mdl_l_t);
      csvwrite([gitRepo '/Data/mdl_l_p.csv'],mdl_l_p);


      [h, crit_p, adj_ci_cvrg, mdl_l_pfdr] = fdr_bh(mdl_l_p,0.001,'pdep','no');
      mdl_l_t(mdl_l_pfdr > 0.01) = 0;
      mdl_l_pfdr(mdl_l_pfdr > 0.01) = 1;

      if subplots == 1
      %plot
      f = figure('units','normalized','outerposition',[0 0 0.55 0.5]);

        ROnSurf = BoSurfStatMakeParcelData(mdl_l_t, G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.0 0.33 0.30 0.30], [0.25 0.33 0.30 0.30], ...
            1, 2, [-max(mdl_l_t) max(mdl_l_t)], cmap_beta);
        hbar = colorbar('horizontal');
        hbar.Position = [0.12 0.1 0.3 0.1];

        ROnSurf = BoSurfStatMakeParcelData((mdl_l_pfdr), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.5 0.33 0.30 0.30], [0.75 0.33 0.30 0.30], ...
            1, 2, [0 1], plasma);
        hbar = colorbar('horizontal');
        hbar.Position = [0.62 0.1 0.3 0.1];

        exportfigbo(f, [gitRepo '/Figures/subplots/AgeChange_PathLenght.png'], 'png', 10)
        close(f)
      end
        %% figure 4
        load([gitRepo '/Utilities/corticalMask.mat'])
        mask = mask.cortex;
        f = figure('units','normalized','outerposition',[0 0 0.8 0.8]);

        ROnSurf = BoSurfStatMakeParcelData(nanmean(Cnorm), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.05 0.7 0.15 0.15], [0.20 0.7 0.15 0.15], ...
            1, 2, [min(nanmean(Cnorm)) max(nanmean(Cnorm))], flipud(plasma));
        hbar = colorbar('horizontal');
        hbar.Position = [0.15 0.7 0.1 0.01];

        ROnSurf = BoSurfStatMakeParcelData(nanmean(Lnorm), G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.4 0.7 0.15 0.15], [0.55 0.7 0.15 0.15], ...
            1, 2, [min(nanmean(Lnorm)) max(nanmean(Lnorm))], flipud(plasma));
        hbar = colorbar('horizontal');
        hbar.Position = [0.50 0.7 0.1 0.01];

        ROnSurf = BoSurfStatMakeParcelData(mdl_c_t, G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.05 0.5 0.15 0.15], [0.20 0.5 0.15 0.15], ...
            1, 2, [-max(mdl_c_t) max(mdl_c_t)], cmap_beta);
        hbar = colorbar('horizontal');
        hbar.Position = [0.15 0.5 0.1 0.01];

        ROnSurf = BoSurfStatMakeParcelData(mdl_l_t, G, parc);
        Rs = SurfStatSmooth(ROnSurf, G, 5);
        Rs = Rs.* mask;
        BoSurfStat_calibrate2Views(Rs, G, [0.4 0.5 0.15 0.15], [0.55 0.5 0.15 0.15], ...
            1, 2, [-max(mdl_l_t) max(mdl_l_t)], cmap_beta);
        hbar = colorbar('horizontal');
        hbar.Position = [0.50 0.5 0.1 0.01];

        exportfigbo(f, [gitRepo '/Figures/Figure_4.png'], 'png', 10)
        close(f)

end
