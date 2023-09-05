function h = xcorr_fft(a,b)
    nans = isnan(a)+isnan(b);
    a(nans>0) = 0;
    b(nans>0) = 0;

    e = fft2(a);
    f = fft2(b);
    
    g = conj(e).*f;
    h = fftshift(ifft2(g));
    
end