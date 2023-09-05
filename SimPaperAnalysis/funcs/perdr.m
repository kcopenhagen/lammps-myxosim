function dxs = perdr(xs, cx, boxSize)
    dxs = xs - cx;
    dxs(dxs>boxSize/2) = dxs(dxs>boxSize/2)-boxSize;
    dxs(dxs<-boxSize/2) = dxs(dxs<-boxSize/2)+boxSize;
end