function adefs = finddefs_sim(bds,boxSize,minr,rsm,XYcal)
    dirf = dfield_sim(bds,boxSize,XYcal,rsm);
    S = sfield_sim(dirf);
    q = chargeeverywhere(dirf);
    
    [adefypix,adefxpix] = find(abs(q)>0.1);
    adefx = adefxpix*XYcal;
    adefy = adefypix*XYcal;
    
    adefind = sub2ind(size(q),adefypix,adefxpix);
    adefq = q(adefind);
    Sdefs = S(adefind);
    
    dx = adefx - adefx';
    dx(dx>boxSize/2) = dx(dx>boxSize/2) - boxSize;
    dx(dx<-boxSize/2) = dx(dx<-boxSize/2) + boxSize;
    dy = adefy - adefy';
    dy(dy>boxSize/2) = dy(dy>boxSize/2) - boxSize;
    dy(dy<-boxSize/2) = dy(dy<-boxSize/2) + boxSize;

    dr = sqrt(dx.^2 + dy.^2);
    
    dr = dr + tril(boxSize*ones(size(dr)));
    
    [mindr, I] = min(dr,[],'all');
    [I,J] = ind2sub(size(dr),I);
    del = [];
    
    while (mindr < minr)
        if (abs(adefq(I) + adefq(J))<0.1)
            del = [del; J; I];
            dr(I,:) = boxSize;
            dr(:,I) = boxSize;
            dr(:,J) = boxSize;
            dr(J,:) = boxSize;
        else
    
            adefq(I) = adefq(I)+adefq(J);
            adefx(I) = (adefx(I)+adefx(J))/2;
            adefy(I) = (adefy(I)+adefy(J))/2;
            del = [del;J];
            dr(:,J) = boxSize;
            adefx(del) = [];
            adefy(del) = [];
            adefq(del) = [];
            dx = adefx - adefx';
            dx(dx>boxSize/2) = dx(dx>boxSize/2) - boxSize;
            dx(dx<-boxSize/2) = dx(dx<-boxSize/2) + boxSize;
            dy = adefy - adefy';
            dy(dy>boxSize/2) = dy(dy>boxSize/2) - boxSize;
            dy(dy<-boxSize/2) = dy(dy<-boxSize/2) + boxSize;
            
            dr = sqrt(dx.^2 + dy.^2);
    
            del = [];
        end

        [mindr, I] = min(dr,[],'all');
        [I,J] = ind2sub(size(dr),I);
    end
    
    adefx(del) = [];
    adefy(del) = [];
    adefq(del) = [];
    
    dirt = dirf;
    dirt(dirt<0) = dirt(dirt<0) + pi;
    dirt = padarray(dirt,[1 1],'circular');
    
    nx = cos(dirt);
    nxnx = nx.*nx;
    ny = sin(dirt);
    nyny = ny.*ny;
    nxny = nx.*ny;
    
    dnxnxdx = diff(nxnx,1,2);
    dnxnydx = diff(nxny,1,2);
    dnxnydy = diff(nxny,1,1);
    dnynydy = diff(nyny,1,1);
    
    dnxnxdx = conv2(dnxnxdx,[0.5; 0.5],'valid');
    dnxnydx = conv2(dnxnydx,[0.5; 0.5],'valid');
    dnxnydy = conv2(dnxnydy,[0.5 0.5],'valid');
    dnynydy = conv2(dnynydy,[0.5 0.5],'valid');
    
    pdefdx = dnxnxdx+dnxnydy;
    pdefdy = dnxnydx+dnynydy;
    pdefdr = sqrt(pdefdx.^2+pdefdy.^2);
    
    pdefdx = pdefdx(2:end,2:end);
    pdefdy = pdefdy(2:end,2:end);
    pdeftheta = atan2(pdefdy,pdefdx);
    
    pdefinds = sub2ind(size(pdeftheta),adefy(adefq>0.1)/XYcal,adefx(adefq>0.1)/XYcal);
    adefd = zeros(size(adefx));
    adefd(adefq>0.1) = pdeftheta(round(pdefinds));
    
    dirt = -dirf;
    dirt(dirt<0) = dirt(dirt<0) + pi;
    dirt = padarray(dirt,[1 1],'circular');
    
    nx = cos(dirt);
    nxnx = nx.*nx;
    ny = sin(dirt);
    nyny = ny.*ny;
    nxny = nx.*ny;
    
    dnxnxdx = diff(nxnx,1,2);
    dnxnydx = diff(nxny,1,2);
    dnxnydy = diff(nxny,1,1);
    dnynydy = diff(nyny,1,1);
    
    dnxnxdx = conv2(dnxnxdx,[0.5; 0.5],'valid');
    dnxnydx = conv2(dnxnydx,[0.5; 0.5],'valid');
    dnxnydy = conv2(dnxnydy,[0.5 0.5],'valid');
    dnynydy = conv2(dnynydy,[0.5 0.5],'valid');
    
    ndefdx = dnxnxdx+dnxnydy;
    ndefdy = dnxnydx+dnynydy;
    ndefdr = sqrt(ndefdx.^2+ndefdy.^2);
    
    ndefdx = ndefdx(2:end,2:end);
    ndefdy = ndefdy(2:end,2:end);
    ndeftheta = atan2(ndefdy,ndefdx);
    
    ndefinds = sub2ind(size(ndeftheta),adefy(adefq<-0.1)/XYcal,adefx(adefq<-0.1)/XYcal);
    adefd(adefq<-0.1) = -ndeftheta(round(ndefinds))/3;
    
    adefs = struct('x',num2cell(adefx),...
        'y',num2cell(adefy),'q',num2cell(adefq),...
        'dx',num2cell(cos(adefd)),'dy',num2cell(sin(adefd)));
end
