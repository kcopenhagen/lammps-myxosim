function l = laserdata(fpath, t)
    %Returns the laser, and height data from path fpath, at time point t. 
    %As well as the time at which the data was taken.
    %Will return NaN if t is beyond the last time point.
    
    name = sprintf('%06d.bin',t-1);

    fname = fullfile(fpath,'Laser',name);
    fid = fopen(fname);
    im = fread(fid,[1024 768],'unsigned short');
    l = transpose(im);
    fclose(fid);

end