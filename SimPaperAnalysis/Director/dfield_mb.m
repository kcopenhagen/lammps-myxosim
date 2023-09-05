function dirf = dfield_mb(img,n)

    % Calculate and return the director field for img smoothing kernel sized n.
        
    % Create gradient kernels for x and y.
    dl = 1/2*[-1 0 1];
    
    dx = zeros(numel(dl),numel(dl));
    dx(end/2+1/2,:) = dl;
    dy = zeros(numel(dl),numel(dl));
    dy(:,end/2+1/2) = dl;
    
    ls = img;
    % Calculate the x and y gradients of the laser image.
    ls = padarray(ls,[1,1],'replicate','both');
    
    lx = conv2(ls,dx,'same');
    ly = conv2(ls,dy,'same');
    
    lx = lx(2:end-1,2:end-1);
    ly = ly(2:end-1,2:end-1);

    % Calculate the hetian for each pixel of the laser image from
    % gradients.
    
    J11 = lx.^2;
    J12 = lx.*ly;
    J22 = ly.^2;

    % Smooth it.
    
    J11 = (imgaussfilt(J11,n));
    J12 = (imgaussfilt(J12,n));
    J22 = (imgaussfilt(J22,n));

    % Calculate the eigenvector for the lower eigenvalue.
    
    VX = -J12;
    VY = (J11-J22)/2+sqrt((J22-J11).^2/4+J12.^2);

    % Save eigenvector direction as the director field.
    
    dirf = atan(VY./VX);
    
end

