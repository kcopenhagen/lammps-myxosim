function [eps, d, l, rho, v, Kagar, Kstiff, rev, run, fp] = getparams(fpath)
   
    if fpath(end) == '/'
        fpath(end) = [];
    end

    folds = strsplit(fpath,'/');
    epss = folds{end-8};
    eps = sscanf(epss,'eps%d');
    ds = folds{end-7};
    d = sscanf(ds,'d%d')/100;
    ls = folds{end-6};
    l = sscanf(ls,'l%d')/100;
    rhos = folds{end-5};
    rho = sscanf(rhos,'rho%d')/100;
    vs = folds{end-4};
    v = sscanf(vs,'v%d')/100;
    Kagars = folds{end-3};
    Kagar = sscanf(Kagars,'Kagar%d');
    Kstiffs = folds{end-2};
    Kstiff = sscanf(Kstiffs,'Kstiff%d');
    revs = folds{end-1};
    rev = sscanf(revs,'rev%d');
    runs = folds{end};
    run = sscanf(runs,'run_%d');

    fp = folds{end-9};
end