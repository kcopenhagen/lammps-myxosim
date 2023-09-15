function dirf = dfield_sim(bds,boxSize,XYcal,rsm)
  %%
  % Note - previously rsm = 7.
  res = ceil(boxSize/XYcal);
  skin = 1;
  ends = zeros(size(bds));
  ends = num2cell(ends);
  [bds.ends] = ends{:};
  %[bds([bds.zs]>0.45).ends] = deal(1);

  %%
  
  xt = [bds([bds.ends]==0).xs]';
  yt = [bds([bds.ends]==0).ys]';
  muxt = [bds([bds.ends]==0).mux]';
  muyt = [bds([bds.ends]==0).muy]';

  xt1 = [xt; ...
      xt(xt>boxSize-skin)-boxSize; xt(xt<skin)+boxSize;...
      xt(yt>boxSize-skin); xt(yt<skin)];
  yt1 = [yt; ...
      yt(xt>boxSize-skin); yt(xt<skin);...
      yt(yt>boxSize-skin)-boxSize; yt(yt<skin)+boxSize];
  muxt1 = [muxt; ...
      muxt(xt>boxSize-skin); muxt(xt<skin);...
      muxt(yt>boxSize-skin); muxt(yt<skin)];
  muyt1 = [muyt; ...
      muyt(xt>boxSize-skin); muyt(xt<skin);...
      muyt(yt>boxSize-skin); muyt(yt<skin)];

  thetas = atan2(muyt1,muxt1);
  thetas(thetas<0) = thetas(thetas<0) + pi;
  muxs = cos(2*thetas);
  muys = sin(2*thetas);

  x = ((boxSize/res)/2):boxSize/res:boxSize;
  y = ((boxSize/res)/2):boxSize/res:boxSize;
  [xx,yy] = meshgrid(x,y);
  
  amux = griddata(xt1,yt1,muxs,xx,yy,'nearest');
  amuy = griddata(xt1,yt1,muys,xx,yy,'nearest');
  
  % amux = medfilt2(amux,[16 16]);
  % amuy = medfilt2(amuy,[16 16]);
  % 
  if (rsm > 0)
    amux = imgaussfilt(amux,rsm);
    amuy = imgaussfilt(amuy,rsm);
  end
  
  dirf = atan2(amuy,amux)/2;

end