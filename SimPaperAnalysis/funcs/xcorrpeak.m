function p = xcorrpeak(h)
    [i,j] = find(h == max(h(:)),1);

    while j<2 || i<2 || j>(numel(h(1,:))-1) || i>(numel(h(:,1))-1)
        h(i,j) = -10;
        [i,j] = find(h == max(h(:)));
    end

    x0 = (log(h(i,j-1)) - log(h(i,j+1)))/(2*log(h(i,j-1))-4*log(h(i,j))...
        +2*log(h(i,j+1)));
    y0 = (log(h(i-1,j)) - log(h(i+1,j)))/(2*log(h(i-1,j))-4*log(h(i,j))...
        +2*log(h(i+1,j)));
    p = [real(j + x0), real(i + y0)];

end