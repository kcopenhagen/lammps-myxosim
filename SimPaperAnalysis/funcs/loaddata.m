function data = loaddata(fpath, t, dataname, datatype)

    if (nargin == 3)
        datatype = dataname;
        dataname = t;
        [fp,fn,ex] = fileparts(fpath);
        fname = fullfile(fp,'..','analysis',dataname,[char(fn) char(ex)]);
        
    else
        name = sprintf('%06d.bin',t-1);
        fname = fullfile(fpath,'analysis',dataname,name);
    end
    fID = fopen(fname,'r');
    data = fread(fID, [768 1024], datatype);
    fclose(fID);
    
end