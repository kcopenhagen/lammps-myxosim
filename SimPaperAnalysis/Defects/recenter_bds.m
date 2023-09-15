function [cx,cy,cmux,cmuy,cvx,cvy] = recenter_bds(bds,def,boxSize,imsz)
    defang = atan2(def.dy,def.dx);
    xs = [bds.xs];
    ys = [bds.ys];
    mux = [bds.mux];
    muy = [bds.muy];
    mur = sqrt(mux.^2 + muy.^2);
    phi = atan2(muy,mux);
    vx = [bds.vx];
    vy = [bds.vy];
    v = sqrt(vx.^2 + vy.^2);
    vang = atan2(vy,vx);
    
    cent_xs = xs - def.x;
    cent_ys = ys - def.y;

    cent_xs(cent_xs<-boxSize/2) = cent_xs(cent_xs<-boxSize/2) + boxSize;
    cent_xs(cent_xs>boxSize/2) = cent_xs(cent_xs>boxSize/2) - boxSize;
    cent_ys(cent_ys<-boxSize/2) = cent_ys(cent_ys<-boxSize/2) + boxSize;
    cent_ys(cent_ys>boxSize/2) = cent_ys(cent_ys>boxSize/2) - boxSize;

    theta = atan2(cent_ys,cent_xs);
    r = sqrt(cent_xs.^2 + cent_ys.^2);
    
    theta = theta - defang;
    phi = phi - defang;
    vang = vang - defang;

    
    cxt = r.*cos(theta);
    cyt = r.*sin(theta);
    cmuxt = mur.*cos(phi);
    cmuyt = mur.*sin(phi);
    cvxt = v.*cos(vang);
    cvyt = v.*sin(vang);
    good = ones(size(cxt));
    good(cxt<-imsz/2) = 0;
    good(cxt>imsz/2) = 0;
    good(cyt<-imsz/2) = 0;
    good(cyt>imsz/2) = 0;
    cx = cxt(good==1);
    cy = cyt(good==1);
    cmux = cmuxt(good==1);
    cmuy = cmuyt(good==1);
    cvx = cvxt(good==1);
    cvy = cvyt(good==1);

end