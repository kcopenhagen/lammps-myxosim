function holes = findholes(l,r,XYcal)
%%
    sz = size(l);
    if numel(sz) == 3
        l = rgb2gray(l);
    end

    l = double(l);
    l = imsharpen(l,'Amount',3,'Radius',3);
    l = imgaussfilt(l,3);
    l = l./imgaussfilt(l,128);
    l = rescale(l,0,1);
    dl = r/(0.66/XYcal) - r/(0.76/XYcal);
    xs = -r/2:r/2;
    ys = -r/2:r/2;
    [xx,yy] = meshgrid(xs,ys);
    rs = sqrt(xx.^2 + yy.^2);
    thering = abs(rs-r/w)<dl;
    
    holes = zeros(size(l));
    
    for xc = 1:dr:(sz(2)-r)
        for yc = 1:dr:(sz(1)-r)
            lt = double(l(yc:yc+r,xc:xc+r));
            lt = abs(fft2(lt));
            lt(1,1) = 0;
        
            lt = fftshift(lt);
            temp = mean(lt(thering));
            holes((yc+r/2-dr/2):(yc+r/2+dr/2),(xc+r/2-dr/2):(xc+r/2+dr/2)) = temp;
        end
    end
    
    holes = double(holes>2);
    holes = round(imgaussfilt(holes,r));
    holes = holes == 0;
    holes = imdilate(holes,strel('disk',dr));
    
end