function [corrl, correrr, fit_corr] = corrlen(drs,corrs,cutoff)
    good = ones(size(drs));
    good(drs<cutoff(1)) = 0;
    good(drs>cutoff(2)) = 0;
    drs_fit = drs(good==1);
    corrs_fit = corrs(good==1);
    fit_corr = fit(drs_fit,corrs_fit,'exp1');    
    errs = confint(fit_corr);

    corrl = -1/fit_corr.b;
    correrr = -1/errs(2,2)+1/fit_corr.b;
end
