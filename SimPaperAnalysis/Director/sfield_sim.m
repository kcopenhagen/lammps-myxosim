function S = sfield_sim(dirf)
    dirf = padarray(dirf,[1 1],'circular');

    nx = cos(2*dirf);
    ny = sin(2*dirf);
    Sx = conv2(nx,0.25*ones(2),'valid');
    Sy = conv2(ny,0.25*ones(2),'valid');

    S = sqrt(Sx.^2 + Sy.^2);
    S = S(2:end,2:end);
end