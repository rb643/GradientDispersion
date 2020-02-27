for s = 1:length(sublist)
    tic
    display(sublist{s});
    lh = load_mgh([inDir '/sub-' sublist{s} '_fmri2fs_bbr_lh_c69-32k.mgh']);
    rh = load_mgh([inDir '/sub-' sublist{s} '_fmri2fs_bbr_lh_c69-32k.mgh']);
    onC69 = [squeeze(lh)' squeeze(rh)'];
    ts_parc = Bo_surfData2parcelData(onC69, SonG);
    R(:,:,s) = corr(ts_parc);
    %Z(:,:,s) = tanh(R);
    dlmwrite([inDir '/func_conn/' sublist{s} '_R.txt'], R(:,:,s))
    toc
end

% if files already exist locally, read them in
parfor s = 1:length(sublist)
   
    R(:,:,s) = dlmread([inDir '/func_conn/' sublist{s} '_R.txt']);
    
end