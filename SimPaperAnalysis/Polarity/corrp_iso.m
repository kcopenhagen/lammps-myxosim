function [rs_iso, corrps_iso] = corrp_iso(fpath,nts,n_pts)

    [files, boxSize] = getframes(fpath);

    res = 1000;
    min_def_dr = 1;

    rs_iso = [];
    corrps_iso = [];

    for n = 1:nts

        ct = randi([300, numel(files)]);

        %%
        bds = loadsimdata2(fullfile(files(ct).folder,files(ct).name));
        dirf = dfield_fast(bds,boxSize,res);
        S = nemorderfield(dirf,12);
        adefs = finddefects(dirf,S,12);
        def_xs = [adefs.x]*boxSize/res;
        def_ys = [adefs.y]*boxSize/res;

        %% Exclude end beads and beads close to defects.
        inds = sub2ind(size(dirf),ceil([bds.ys]*res/boxSize),ceil([bds.xs]*res/boxSize));
        ends = zeros(size(bds));
        muxs = [bds.mux];
        muys = [bds.muy];
        mumags = sqrt(muxs.^2 + muys.^2);
        muxs = muxs./mumags;
        muys = muys./mumags;

        ends(abs(muxs.*cos(dirf(inds)) + muys.*sin(dirf(inds))) < 0.96) = 1;

        defdxs = [bds.xs] - def_xs';
        defdys = [bds.ys] - def_ys';
        defdrs = sqrt(defdxs.^2 + defdys.^2);
        mindefdrs = min(defdrs);
        ends(mindefdrs<min_def_dr) = 1;

        ends = num2cell(ends);
        [bds.ends] = ends{:};

        bds([bds.ends]==1) = [];
        %%

        inds = randperm(numel(bds));

        ind1s = inds(1:n_pts);
        ind2s = inds(n_pts+1:2*n_pts);

        dxs = [bds(ind1s).xs] - [bds(ind2s).xs]';
        dys = [bds(ind1s).ys] - [bds(ind2s).ys]';
        drs = sqrt(dxs.^2 + dys.^2);

        muxs1 = [bds(ind1s).mux];
        muys1 = [bds(ind1s).muy];

        muxs2 = [bds(ind2s).mux];
        muys2 = [bds(ind2s).muy];

        pcorrs = muxs1.*muxs2' + muys1.*muys2';

        rs_iso = [rs_iso; drs(:)];
        corrps_iso = [corrps_iso; pcorrs(:)];
        save('corrpiso_inprog.mat','rs_iso','corrps_iso')
    end

end