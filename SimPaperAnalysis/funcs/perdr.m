function dxs = perdr(xs1, xs2, boxSize)
    dxs = xs1 - xs2';
    dxs(dxs>boxSize/2) = dxs(dxs>boxSize/2)-boxSize;
    dxs(dxs<-boxSize/2) = dxs(dxs<-boxSize/2)+boxSize;
end