function adefs = findefs_sim(bds,r,XYcal,boxSize)
    res = ceil(boxSize/XYcal);
    resz = round(r/XYcal);
    [re, rimsz] = makeringelement(resz);
    dirf = dfield_sim(bds,boxSize,res,0);
    
    qx = 1:res;
    qy = 1:res;
    [qxx, qyy] = meshgrid(qx,qy);
    qxx = qxx(:);
    qyy = qyy(:);
    aqs = zeros(size());
    for i = 1:numel(qxx)
        dirt = dirf((qxx(i)-(rimsz-1)/2):(qxx(i)+(rimsz-1)/2),...
            (qyy(i)-(rimsz-1)/2):(qyy(i)+(rimsz-1)/2));
        angs = dirt(re);
        dangs = diff(angs);
        dangs(dangs>pi/2) = dangs(dangs>pi/2)-pi;
        dangs(dangs<-pi/2) = dangs(dangs<-pi/2)+pi;
        q = sum(dangs);
        
    end
    
end
