function data = loaddata(fpath, t, dataname, datatype)

    name = sprintf('%06d.bin',t-1);
    fID = fopen([fpath '/analysis/' dataname '/' name],'r');
    data = fread(fID, [768 1024], datatype);
    fclose(fID);
    
end