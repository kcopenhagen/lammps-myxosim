function polf = pfield_sim(bds,boxSize,res)

  skin = 1;

  ends = zeros(size(bds));
  ends = num2cell(ends);
  [bds.ends] = ends{:};
  [bds([bds.zs]>0.45).ends] = deal(1);

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
  %thetas(thetas<0) = thetas(thetas<0) + pi;
  muxs = cos(thetas);
  muys = sin(thetas);

  x = boxSize/res:boxSize/res:boxSize;
  y = boxSize/res:boxSize/res:boxSize;
  [xx,yy] = meshgrid(x,y);
  amux = griddata(xt1,yt1,muxs,xx,yy,'nearest');
  amuy = griddata(xt1,yt1,muys,xx,yy,'nearest');
  
  polf = atan2(amuy,amux);

end