function laservid(fpath,stT,finT)
  res = 1389;
  addpath(genpath("~/SimPaperAnalysis"));
  
  [files,boxSize] = getframes(fpath);

  v = VideoWriter(fullfile(fpath,'laservid.avi'),'Motion JPEG AVI');
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
    
    cellimg = cellim(bds,boxSize,res,[0 0],4,2,1);
    cellimg(cellimg>0.7) = 0.7;
    rescale(cellimg,0,1);
    cellimg = real2rgb(cellimg,gray);
    writeVideo(v,cellimg);
  end
  close(v)
  disp("laser vid complete")
end
