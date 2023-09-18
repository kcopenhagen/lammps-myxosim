
function [cx,cy,cvx,cvy] = recenter_eVs(Vx,Vy,def,XYcal,imsz)

    defx = def.x/XYcal;
    defy = def.y/XYcal;
    defang = atan2(def.dy,def.dx);

    sz = size(Vx);
    xxs = (1:sz(2))*XYcal;
    yys = (1:sz(1))*XYcal;
    [xs,ys] = meshgrid(xxs,yys);
    xs = xs(:);
    ys = ys(:);
    Vx = Vx(:);
    Vy = Vy(:);

    v = sqrt(Vx.^2 + Vy.^2);
    vang = atan2(Vy,Vx);
    
    cent_xs = xs - def.x;
    cent_ys = ys - def.y;

    theta = atan2(cent_ys,cent_xs);
    r = sqrt(cent_xs.^2 + cent_ys.^2);
    
    theta = theta - defang;
    vang = vang - defang;
    
    cxt = r.*cos(theta);
    cyt = r.*sin(theta);
    cvxt = v.*cos(vang);
    cvyt = v.*sin(vang);
    good = ones(size(cxt));
    good(cxt<-imsz/2) = 0;
    good(cxt>imsz/2) = 0;
    good(cyt<-imsz/2) = 0;
    good(cyt>imsz/2) = 0;
    cx = cxt(good==1);
    cy = cyt(good==1);
    cvx = cvxt(good==1);
    cvy = cvyt(good==1);
end