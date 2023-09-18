function [rs_par, rs_perp, corrp_par, corrp_perp] = corrp_parperp(fpath, n_ts, n_pts)
%%
    [files, boxSize] = getframes(fpath);
    res = 1000;
    dr = 0.05;
    max_dr = 5;
    max_L = boxSize*0.8;
    max_cl = 0.24;

    min_def_dr = 1;

    %% 

    rs_par = [];
    rs_perp = [];
    corrp_par = [];
    corrp_perp = [];

    %%
    for n = 1:n_ts
        %%
        ct = randi([1, numel(files)]);
        %%
        bds = loadsimdata2(fullfile(files(ct).folder,files(ct).name));
        dirf = dfield_fast(bds,boxSize,res);
        S = nemorderfield(dirf,12);
        adefs = finddefects(dirf,S,12);
        def_xs = [adefs.x]*boxSize/res;
        def_ys = [adefs.y]*boxSize/res;
%%
        inds = sub2ind(size(dirf),ceil([bds.ys]*res/boxSize),ceil([bds.xs]*res/boxSize));
        ends = zeros(size(bds));
        muxs = [bds.mux];
        muys = [bds.muy];
        mumags = sqrt(muxs.^2 + muys.^2);
        muxs = muxs./mumags;
        muys = muys./mumags;

        ends(abs(muxs.*cos(dirf(inds)) + muys.*sin(dirf(inds))) < 0.96) = 1;
        ends = num2cell(ends);
        [bds.ends] = ends{:};

        %%
        for i = 1:n_pts
          %%
          pt = randi(numel(bds));
          
          % Regenerate if the angle between the dfield and polarity vector is large.
          % (Should only be big for lead beads which have noise).
          while (bds(pt).ends == 1)
              pt = randi(numel(bds));
             
          end
          cx = bds(pt).xs;
          cy = bds(pt).ys;
          cmux = cos(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
              ceil(mod(cx,boxSize)*res/boxSize)));
          cmuy = sin(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
              ceil(mod(cx,boxSize)*res/boxSize)));
          
          go = 1;
          pmux = cmux;
          pmuy = cmuy;
          [smux,smuy] = unit_v(bds(pt).mux,bds(pt).muy);
          cpt = pt;
          tdr = 0;
          rs_par = [rs_par; tdr];
          corrp_par =[corrp_par; smux*smux + smuy*smuy];
          
          sx = cx;
          sy = cy;
          bds = recentbds(bds,cx,cy,boxSize);
          def_xs((def_xs-cx) > boxSize/2) = def_xs((def_xs-cx) > boxSize/2) - boxSize;
          def_xs((def_xs-cx) <-boxSize/2) = def_xs((def_xs-cx) <-boxSize/2) + boxSize;

          def_ys((def_ys-cy) > boxSize/2) = def_ys((def_ys-cy) > boxSize/2) - boxSize;
          def_ys((def_ys-cy) <-boxSize/2) = def_ys((def_ys-cy) <-boxSize/2) + boxSize;

          fwds = 1;
          %%
          % Calculate p correlation along director field.
          npts = pt;
          sl_x = [];
          sl_y = [];
          while (go == 1)
              if pmux*cmux + pmuy*cmuy < 0
                  cmux = -cmux;
                  cmuy = -cmuy;
              end

              nx = cx + dr * cmux;
              ny = cy + dr * cmuy;
              def_dx = def_xs - nx;
              def_dy = def_ys - ny;
              def_drs = sqrt(def_dx.^2 + def_dy.^2);
              if (min(def_drs) < min_def_dr)
                  go = 0;
              end
              sl_x = [sl_x; nx];
              sl_y = [sl_y; ny];
              tdr = tdr + dr;

              dxs = nx - [bds.xs];
              dys = ny - [bds.ys];
              drs = sqrt(dxs.^2 + dys.^2);
              [mind, npt] = min(drs);

              cx = nx;
              cy = ny;
              pmux = cmux;
              pmuy = cmuy;
              cmux = cos(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                  ceil(mod(cx,boxSize)*res/boxSize)));
              cmuy = sin(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                  ceil(mod(cx,boxSize)*res/boxSize)));
              if (mind < max_cl && bds(npt).ends == 0 && npt ~= cpt && go == 1)
                  tdr = tdr + mind;
                  cx = bds(npt).xs;
                  cy = bds(npt).ys;
                  cmux = cos(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                      ceil(mod(cx,boxSize)*res/boxSize)));
                  cmuy = sin(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                      ceil(mod(cx,boxSize)*res/boxSize)));
                  
                  cpt = npt;

                  rs_par = [rs_par; tdr];
                  [cmuxp, cmuyp] = unit_v(bds(cpt).mux,bds(cpt).muy);
                  corrp_par = [corrp_par; smux*cmuxp + smuy*cmuyp];
                  npts = [npts; cpt];
              end
              if (tdr > max_L/2)
                  go = 0;
              end
              if (tdr - rs_par(end) > max_dr)
                  go = 0;
              end
              if (go == 0 && fwds == 1)
                  fwds = 0;
                  cx = sx;
                  cy = sy;
                  cmux = cos(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                      ceil(mod(cx,boxSize)*res/boxSize)));
                  cmuy = sin(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                      ceil(mod(cx,boxSize)*res/boxSize)));
                  pmux = -cmux;
                  pmuy = -cmuy;
                  go = 1;
                  tdr = 0;
                  cpt = pt;
              end
          end
          %% Pcorr along perpendicular direction

          cx = bds(pt).xs;
          cy = bds(pt).ys;

          cmux = cos(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
              ceil(mod(cx,boxSize)*res/boxSize))+pi/2);
          cmuy = sin(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
              ceil(mod(cx,boxSize)*res/boxSize))+pi/2);
          
          go = 1;
          pmux = cmux;
          pmuy = cmuy;
          [smux, smuy] = unit_v(bds(pt).mux,bds(pt).muy);
          cpt = pt;
          tdr = 0;
          rs_perp = [rs_perp; tdr];
          corrp_perp =[corrp_perp; smux*smux + smuy*smuy];
          
          sx = cx;
          sy = cy;
          bds = recentbds(bds,cx,cy,boxSize);
          def_xs((def_xs-cx) > boxSize/2) = def_xs((def_xs-cx) > boxSize/2) - boxSize;
          def_xs((def_xs-cx) <-boxSize/2) = def_xs((def_xs-cx) <-boxSize/2) + boxSize;

          def_ys((def_ys-cy) > boxSize/2) = def_ys((def_ys-cy) > boxSize/2) - boxSize;
          def_ys((def_ys-cy) <-boxSize/2) = def_ys((def_ys-cy) <-boxSize/2) + boxSize;

          fwds = 1;

          nptsp = pt;
          sl_xp = [];
          sl_yp = [];
          while (go == 1)
              if pmux*cmux + pmuy*cmuy < 0
                  cmux = -cmux;
                  cmuy = -cmuy;
              end

              nx = cx + dr * cmux;
              ny = cy + dr * cmuy;
              def_dx = def_xs - nx;
              def_dy = def_ys - ny;
              def_drs = sqrt(def_dx.^2 + def_dy.^2);
              if (min(def_drs) < min_def_dr)
                  go = 0;
              end
              sl_xp = [sl_xp; nx];
              sl_yp = [sl_yp; ny];
              tdr = tdr + dr;

              dxs = nx - [bds.xs];
              dys = ny - [bds.ys];
              drs = sqrt(dxs.^2 + dys.^2);
              [mind, npt] = min(drs);

              cx = nx;
              cy = ny;
              pmux = cmux;
              pmuy = cmuy;
              cmux = cos(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                  ceil(mod(cx,boxSize)*res/boxSize))+pi/2);
              cmuy = sin(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                  ceil(mod(cx,boxSize)*res/boxSize))+pi/2);
              if (mind < max_cl && bds(npt).ends == 0 && npt ~= cpt && go == 1)
                  tdr = tdr + mind;
                  cx = bds(npt).xs;
                  cy = bds(npt).ys;
                  cmux = cos(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                      ceil(mod(cx,boxSize)*res/boxSize))+pi/2);
                  cmuy = sin(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                      ceil(mod(cx,boxSize)*res/boxSize))+pi/2);
                  
                  cpt = npt;

                  [cmuxp, cmuyp] = unit_v(bds(cpt).mux,bds(cpt).muy);
                  rs_perp = [rs_perp; tdr];
                  corrp_perp = [corrp_perp; smux*cmuxp + smuy*cmuyp];
                  nptsp = [nptsp; cpt];
              end
              if (tdr > max_L/2)
                  go = 0;
              end
              if (tdr - rs_par(end) > max_dr)
                  go = 0;
              end
              if (go == 0 && fwds == 1)
                  fwds = 0;
                  cx = sx;
                  cy = sy;
                  cmux = cos(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                      ceil(mod(cx,boxSize)*res/boxSize))+pi/2);
                  cmuy = sin(dirf(ceil(mod(cy,boxSize)*res/boxSize),...
                      ceil(mod(cx,boxSize)*res/boxSize))+pi/2);
                  pmux = -cmux;
                  pmuy = -cmuy;
                  go = 1;
                  tdr = 0;
                  cpt = pt;
              end
          end
	    save('corrpparperp_inprog.mat','rs_par','rs_perp','corrp_par','corrp_perp')
        end
    end
end

