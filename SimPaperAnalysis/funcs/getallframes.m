function files = getallframes(dataset)
switch dataset
    case "sim"
        eps = 5;
        d = 0.7;
        l = 7;
        rho = 0.3;
        v = 5;
        Kagar = 500;
        Kstiff = 100;
        rev = 8;
        runs = 0:2;
        fp = "PoissonRevs";
        files = [];
    
        for r = 1:numel(runs)
            run = runs(r);
            fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,run,fp);
            [filest, ~] = getframes(fpath);
            for f = 1:numel(filest)
                if filest(f).name(1) ~= '.'
                    file = fullfile(filest(f).folder,filest(f).name);
                    files = [files; string(file)];
                end
            end
        end
    case "mbwt"
        [fpath1,fpath2,fpath3] = exptpaths("mbwt");
        files = [];
        for f = 1:numel(fpath3)
            fpath = fullfile(fpath1,fpath2,fpath3(f),'img');
            filest = dir(fpath);
            filest = filest(~[filest.isdir]);
            for g = 1:numel(filest)
                if filest(g).name(1) ~= '.'
                    files = [files; string(fullfile(filest(g).folder,...
                        filest(g).name))];
                end
            end
        end
        
    case "mbpila"
        [fpath1,fpath2,fpath3] = exptpaths("mbpila");

        files = [];
        for f = 1:numel(fpath3)
            fpath = fullfile(fpath1,fpath2,fpath3(f),'img');
            filest = dir(fpath);
            filest = filest(~[filest.isdir]);
            for g = 1:numel(filest)
                if filest(g).name(1) ~= '.'
                    files = [files; string(fullfile(filest(g).folder,...
                        filest(g).name))];
                end
            end
        end
        
    case "kcwt"
        [fpath1,fpath2,~] = exptpaths("kcwt");

        files = [];
        for f = 1:numel(fpath2)
            fpath = fullfile(fpath1,fpath2(f),'Laser');
            filest = dir(fpath);
            filest = filest(~[filest.isdir]);
            for g = 1:numel(filest)
                if filest(g).name(1) ~= '.'
                    files = [files; string(fullfile(filest(g).folder,...
                        filest(g).name))];
                end
            end
        end

    otherwise
        disp("Error getting file paths");

end



