function holes = findholes(img,n)

    % Calculate and return the director field for img smoothing kernel sized n.
        
    % Preprocess img.
    sz = size(img);
    if numel(sz)==3
        img = rgb2gray(img);
    end
    img = double(img);
    img = imsharpen(img,'Amount',3,'Radius',3);
    img = imgaussfilt(img,3);
    img = img./imgaussfilt(img,128);
    img = rescale(img,0,1);
    img = histeq(img,200);

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

    % Eigenvalues.
    l1 = ((J11+J22) + sqrt((J11+J22).^2 - 4*(J11.*J22-J12.^2)))/2;
    l2 = ((J11+J22) - sqrt((J11+J22).^2 - 4*(J11.*J22-J12.^2)))/2;
    
    ldiff = l1-l2;
    holes = ldiff<0.002;
    holes = imopen(holes,strel('disk',n+4));
    holes = imdilate(holes,strel('disk',n));
end
