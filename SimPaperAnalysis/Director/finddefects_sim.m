function adefs = finddefects_sim(dfield, S, r, boxSize, res)
%Finds and returns the x, and y coordinates of defects the charge and the
%direction of each defect at time t, in experiment fpath, as arrays.

    % Load in director field and order field for the current frame. Find
    % the ring element for the radius r. (circle of pixels).
    % r = 6;
    re = makeringelement(r);
    %S = loaddata(fpath,t,'order','float');
    %dfield = loaddata(fpath,t,'dfield','float');
    
    defs = (S<0.5);
    %defs = imextendedmin(S,0.3).*imextendedmin(S,0.2);
    
    CC = bwconncomp(defs);
    P = regionprops(CC, 'Centroid');
    Centroids = [P.Centroid];
    Centx = round(Centroids(1:2:end-1));
    Centy = round(Centroids(2:2:end));
    %Exclude defects too close to the wall to define the ring to calculate
    %charge with.
    nog = ~(Centx<r).*~(Centx>numel(S(1,:))-r).*~(Centy<r).*~(Centy>numel(S(:,1))-r);
    Centx(~nog) = [];
    Centy(~nog) = [];
    
    q = [];
    d = [];
    
    %Loop through all defects to find charge and direction for them.
    for i = 1:numel(Centx)
        %Get the director field just around the defect centroid.
        dirt = dfield(Centy(i)-r+1:Centy(i)+r-1,Centx(i)-r+1:Centx(i)+r-1);
        %Find the director at each ring element pixel.
        ang1s = dirt(re);
        ang2s = [ang1s(2:end); ang1s(1)];
        ang1s(isnan(ang1s))=[];
        ang2s(isnan(ang2s))=[];
        %Calculate the change in angle from pixel to pixel. Correct for
        %looping around.
        das = ang2s-ang1s;
        das(das>pi/2)=das(das>pi/2)-pi/2;
        das(das<-pi/2)=das(das<-pi/2)+pi/2;
        %Set charge of defect to sum of angle changes over pi.
        q = [q sum(das)/pi];
        %Calculate direction of positive defects.
        if abs(q(i)-0.5) < 0.1
            nx = cos(dirt);
            nxnx = nx.*nx;
            ny = sin(dirt);
            nyny = ny.*ny;
            nxny = nx.*ny;
            [dnxnxdx, ~] = gradient(nxnx);
            [dnxnydx, dnxnydy] = gradient(nxny);
            [~, dnynydy] = gradient(nyny);
            dnxnxdx = dnxnxdx(r,r);
            dnxnydx = dnxnydx(r,r);
            dnxnydy = dnxnydy(r,r);
            dnynydy = dnynydy(r,r);
            dx = dnxnxdx+dnxnydy;
            dy = dnxnydx+dnynydy;
            dr = sqrt(dx^2+dy^2);
            d = [d; [dx/dr dy/dr]];
            
        %Calculate direction of negative defects.
        elseif abs(q(i)+0.5) < 0.1
            
            nx = cos(-dirt);
            nxnx = nx.*nx;
            ny = sin(-dirt);
            nyny = ny.*ny;
            nxny = nx.*ny;
            
            [dnxnxdx, ~] = gradient(nxnx);
            [dnxnydx, dnxnydy] = gradient(nxny);
            [~, dnynydy] = gradient(nyny);
            dnxnxdx = dnxnxdx(r,r);
            dnxnydx = dnxnydx(r,r);
            dnxnydy = dnxnydy(r,r);
            dnynydy = dnynydy(r,r);
            
            psiprime = atan2(dnxnydx+dnynydy,dnxnxdx+dnxnydy);
            d = [d; [cos(-psiprime/3) sin(-psiprime/3)]];
        else
            d = [d; NaN NaN];
        end
    end
    
    x = Centx*boxSize/res;
    y = Centy*boxSize/res;
    
    if numel(x)>0
        adefs = struct('x',num2cell(x'),'y',num2cell(y'),'q',num2cell(q'),'d',num2cell(d,2));
        %adefs(del) = [];
        adefs(abs([adefs.q])<0.1) = [];
    else
        adefs = [];
    end
    
end
    

function inds = makeringelement(r)
    seb = strel('disk',r);
    sea = strel('disk',r-1);
    sea = padarray(sea.Neighborhood,[1 1]);
    seneigh = seb.Neighborhood - sea;
    [y, x] = find(seneigh);
    x = x-numel(seneigh(:,1))/2;
    y = y-numel(seneigh(:,1))/2;
    a = atan2(y,x);
    [~,ind] = sort(a);
    x = x(ind);
    y = y(ind);

    x = x+numel(seneigh(:,1))/2;
    y = y+numel(seneigh(:,1))/2;
    inds = sub2ind(size(seneigh),y,x);
end
    