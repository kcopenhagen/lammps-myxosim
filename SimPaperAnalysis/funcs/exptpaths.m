
function [fpath1, fpath2, fpath3] = exptpaths(ename)
%% Returns paths to different experiments:
% mbpila = matt's pilA data.
% mbwt = matt's wt data.
% kcwt = katie's wt data.

if nargin < 1
    disp("acceptable inputs - mbpila, mbwt, kcwt.")
else
%% Matt's data.
    if ename == "mbpila"
        fpath1 = '/tigress/SHAEVITZ/mb46/rgb_cam/data/myxo/';
        
        %% pilA
        
        fpath2 = 'dk_pilA';
        
        fpath3 = ["2022-06-16/trial2"; 
            "2022-06-16/trial2c";
            "2022-06-16/trial3";
            "2022-06-23/trial2";
            "2022-06-23/trial3";
            "2022-06-23/trial4";
            "2022-10-31/trial1";
            "2022-10-31/trial3";
            "2022-11-08/trial1";
            "2022-11-08/trial2";];
    %% WT
    elseif ename=="mbwt"
        fpath1 = '/tigress/SHAEVITZ/mb46/rgb_cam/data/myxo/';
        fpath2 = 'dk_wt';
        fpath3 = ["2022-06-18/trial1";
            "2022-06-18/trial2";
            "2022-06-18/trial3";
            "2022-06-21/trial2";];
    %% My data
    elseif ename == "kcwt"
        fpath1 = '/tigress/kc32/Monolayerpaperdata/Data/';
        fpath2 = ["190109KC3";
            "190109KC3b";
            "190111KC2";
            "190116KC2";
            "190117KC1";
            "190117KC2";
            "190117KC2a";
            "190117KC2b";];
        fpath3 = [];
    else
        disp("acceptable inputs - mbpila, mbwt, kcwt.")
    end
end

