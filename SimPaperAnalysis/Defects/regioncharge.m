function q = regioncharge(dirf,Reg,method)
    
    Reg2 = imdilate(Reg,strel('disk',1));
    outline = Reg2-Reg;

    [Ys, Xs] = find(outline==1);
    COMx = mean(Xs);
    COMy = mean(Ys);
    angs = atan2(Ys-COMy,Xs-COMx);

    if method == "fast"
        [~,inds] = sort(angs);
        Xs = Xs(inds);
        Ys = Ys(inds);
        sinds = sub2ind(size(Reg),Ys,Xs);
        sinds = [sinds; sinds(1)];
        dirfR = dirf(sinds);
        diffs = diff(dirfR);
        diffs(diffs>pi/2) = diffs(diffs>pi/2) - pi;
        diffs(diffs<-pi/2) = diffs(diffs<-pi/2) + pi;
        q = sum(diffs)/(2*pi);
    elseif method == "slow"
        dxs = Xs - Xs';
        dys = Ys - Ys';
        drs = dxs.^2 + dys.^2;
        inds = 1;
        drs(:,1) = 999;
        while (numel(inds)<numel(Xs))
            drs(drs<0.5) = 999;
            [~,j] = min(drs(inds(end),:));
            drs(:,j) = 999;
            inds = [inds; j];
        end
        Xs = Xs(inds);
        Ys = Ys(inds);
        sangs = angs(inds);
        
        sinds = sub2ind(size(Reg),Ys,Xs);
        sinds = [sinds; sinds(1)];
        dirfR = dirf(sinds);
        
        diffsangs = diff(sangs);
        diffsangs(diffsangs>pi/2) = diffsangs(diffsangs>pi/2) - pi;
        diffsangs(diffsangs<-pi/2) = diffsangs(diffsangs<-pi/2) + pi;

        diffs = diff(dirfR);
        diffs(diffs>pi/2) = diffs(diffs>pi/2) - pi;
        diffs(diffs<-pi/2) = diffs(diffs<-pi/2) + pi;
        q = sum(diffs)/sum(diffsangs);

    end

end
    