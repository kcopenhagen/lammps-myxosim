function data = loaddata(fpath, t, dataname, datatype)

    name = sprintf('%06d.bin',t-1);
    fname = fullfile(fpath,'analysis',dataname,name);
    fID = fopen(fname,'r');
    data = fread(fID, [768 1024], datatype);
    fclose(fID);
    
end