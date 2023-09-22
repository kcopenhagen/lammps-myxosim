function adefs = finddefs_sim(bds,boxSize,minr,rsm,XYcal,dsm)
    % 

    dirf = dfield_sim(bds,boxSize,XYcal,rsm);
    %S = sfield_sim(dirf);
    q = chargeeverywhere(dirf);
    
    [adefypix,adefxpix] = find(abs(q)>0.1);
    adefx = adefxpix*XYcal;
    adefy = adefypix*XYcal;
    
    adefind = sub2ind(size(q),adefypix,adefxpix);
    adefq = q(adefind);
    %Sdefs = S(adefind);
    
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
            del = [del; J];
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
            dr = dr + tril(boxSize*ones(size(dr)));

            del = [];
        end

        [mindr, I] = min(dr,[],'all');
        [I,J] = ind2sub(size(dr),I);
    end
    
    adefx(del) = [];
    adefy(del) = [];
    adefq(del) = [];
    
    dirt = dirf+pi/4;
    dirt = padarray(dirt,[1 1],'circular');
    
    nx = cos(dirt);
    ny = sin(dirt);
    
    nxnx = nx.*nx;
    nxny = nx.*ny;
    nyny = ny.*ny;
    
    dxfilt = [-1 0; 0 1];
    dyfilt = [0 -1; 1 0];
    dnxnxdx = conv2(nxnx,dxfilt,'valid');
    dnxnydx = conv2(nxny,dxfilt,'valid');
    dnxnydy = conv2(nxny,dyfilt,'valid');
    dnynydy = conv2(nyny,dyfilt,'valid');
    
    if (dsm>0)
        dsmpix = dsm/XYcal;
        dnxnxdx = imgaussfilt(dnxnxdx,dsmpix,'Padding','circular');
        dnxnydx = imgaussfilt(dnxnydx,dsmpix,'Padding','circular');
        dnxnydy = imgaussfilt(dnxnydy,dsmpix,'Padding','circular');
        dnynydy = imgaussfilt(dnynydy,dsmpix,'Padding','circular');
    end

    pdefdx = dnxnxdx+dnxnydy;
    pdefdy = dnxnydx+dnynydy;
    
    pdefdx = pdefdx(2:end,2:end);
    pdefdy = pdefdy(2:end,2:end);
    pdeftheta = atan2(pdefdy,pdefdx)+pi/4;

    pdefinds = sub2ind(size(pdeftheta),adefy(adefq>0.1)/XYcal,adefx(adefq>0.1)/XYcal);
    adefd = zeros(size(adefx));
    adefd(adefq>0.1) = pdeftheta(round(pdefinds));
    
    dirt = -(dirf+pi/4);

    dirt = padarray(dirt,[1 1],'circular');
    
    nx = cos(dirt);
    ny = sin(dirt);
    
    nxnx = nx.*nx;
    nxny = nx.*ny;
    nyny = ny.*ny;

    dnxnxdx = conv2(nxnx,dxfilt,'valid');
    dnxnydx = conv2(nxny,dxfilt,'valid');
    dnxnydy = conv2(nxny,dyfilt,'valid');
    dnynydy = conv2(nyny,dyfilt,'valid');

    if (dsm>0)
        dsmpix = dsm/XYcal;
        dnxnxdx = imgaussfilt(dnxnxdx,dsmpix,'Padding','circular');
        dnxnydx = imgaussfilt(dnxnydx,dsmpix,'Padding','circular');
        dnxnydy = imgaussfilt(dnxnydy,dsmpix,'Padding','circular');
        dnynydy = imgaussfilt(dnynydy,dsmpix,'Padding','circular');
    end

    ndefdx = dnxnxdx+dnxnydy;
    ndefdy = dnxnydx+dnynydy;
    
    ndefdx = ndefdx(2:end,2:end);
    ndefdy = ndefdy(2:end,2:end);
    ndeftheta = atan2(ndefdy,ndefdx)+pi/4;
    
    ndefinds = sub2ind(size(ndeftheta),adefy(adefq<-0.1)/XYcal,adefx(adefq<-0.1)/XYcal);
    adefd(adefq<-0.1) = -ndeftheta(round(ndefinds))/3+pi/3;
    
    adefs = struct('x',num2cell(adefx),...
        'y',num2cell(adefy),'q',num2cell(adefq),...
        'dx',num2cell(cos(adefd)),'dy',num2cell(sin(adefd)));

end
