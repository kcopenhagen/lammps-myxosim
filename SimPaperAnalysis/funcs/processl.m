function l = processl(lraw,dataset)    
    sz = size(lraw);
    if (numel(sz) == 3)
        lt = rgb2gray(lraw);
    else
        lt = lraw;
    end
    lt = double(imsharpen(lt,'Amount',3,'Radius',3));
    switch dataset
        case {"mbwt", "mbpila", "mbflavo"}
            lt = imgaussfilt(lt,3);
        case "kcwt"
            lt = imgaussfilt(lt,1);
        otherwise
            lt = lt;
    end
    switch dataset
        case "sim"
            l = lraw;
        otherwise
            lt = rescale(lt,0,1);
            % lt = lt./imgaussfilt(lt,128);
            lt = adapthisteq(lt);
            l = adapthisteq(lt);
    end
end