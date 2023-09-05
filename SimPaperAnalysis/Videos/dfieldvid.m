function dfieldvid(fpath,stT,finT)
  res = 1000;
  addpath(genpath("~/SimPaperAnalysis"));
  
  [files,boxSize] = getframes(fpath);

  v = VideoWriter(fullfile(fpath,'dfieldvid.avi'),'Motion JPEG AVI');
  v.Quality = 95;
  v.FrameRate = 12;
  open(v);
  stT = 1;
  finT = numel(files);

  for t = stT:finT
    if (mod(t,10)==0)
        disp(sprintf('t = %d',t));
    end
    fname = fullfile(files(t).folder,files(t).name);
    bds = loadsimdata(fname);
    dirf = dfield_sim(bds,boxSize,res);
    img = real2rgb(dirf,orientcmap);
%%
    sz2 = size(dirf);
    cellimg = zeros(sz2);

    indx = round([bds.xs]*res/boxSize)+1;
    indy = round([bds.ys]*res/boxSize)+1;

    gd = ones(size(indx));
    gd(indx>sz2(2)) = 0;
    gd(indy>sz2(1)) = 0;

    indx = indx(gd==1);
    indy = indy(gd==1);

    inds = sub2ind(sz2,indy,indx);
    cellimg(inds) = 1;
    cellimg = imdilate(cellimg, strel('disk',round(0.3*res/boxSize)));
    cellimg = imgaussfilt(cellimg, round(0.2*res/boxSize));
    cellimg = histeq(cellimg);

    imgr = img(:,:,1);
    imgg = img(:,:,2);
    imgb = img(:,:,3);

    imgr = imgr.*(1-0.25*(1-cellimg));
    imgg = imgg.*(1-0.25*(1-cellimg));
    imgb = imgb.*(1-0.25*(1-cellimg));

    img(:,:,1) = imgr;
    img(:,:,2) = imgg;
    img(:,:,3) = imgb;

    writeVideo(v,img);

  end
  close(v)
  disp("pfield vid complete")
end
