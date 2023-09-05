function img = loadrgbimg(fname,norm)
    img = importdata(fname);
    img = im2double(rgb2gray(img));
    if nargin==1
        norm = false;
    end
    if norm
        img = img./imgaussfilt(img,128);
        img = imgaussfilt(img,1);
    end
end