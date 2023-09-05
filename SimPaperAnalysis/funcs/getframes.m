
function [files, boxSize] = getframes(fpath)

    files = dir(fullfile(fpath,"cells"));
    files(arrayfun(@(x)numel(x.name)<3,files)) = [];
    datafiles = arrayfun(@(x) x.name(end-2:end)=="txt",files);
    files = files(datafiles);
    fid = fopen(fullfile(fpath,'spfr_init.txt'));
    boxSize = [];
    while isempty(boxSize)
        tline = fgetl(fid);
        if contains(tline,'xlo')
            tt = strsplit(tline);
            boxSize = str2num(tt{3}) - str2num(tt{2});
        end
    end
    fclose(fid);
end