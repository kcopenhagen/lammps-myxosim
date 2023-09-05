function ts = getts(fpath)
    fid = fopen(fullfile(fpath,'times.txt'));
    ts = fscanf(fid,'%f');
    fclose(fid);
    
    files = dir(fullfile(fpath,'Laser/'));
    dF = [files.bytes]>1500000;
    files = files(dF);
    ts = ts(1:numel(files));

end