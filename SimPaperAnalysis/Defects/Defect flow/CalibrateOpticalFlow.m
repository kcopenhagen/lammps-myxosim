
dataset = "mbwt";
disp(dataset);

switch dataset
    case "mbwt"
        [fpath1, fpath2, fpath3] = exptpaths(dataset);
    case "mbpila"
        [fpath1, fpath2, fpath3] = exptpaths(dataset);
    case "kcwt"
        [fpath1, fpath2, ~] = exptpaths(dataset);
    case "sim"
        
    otherwise
end

opticFlow = opticalFlowFarneback('NumPyramidLevels',3,'NumIterations',1,...
    'NeighborhoodSize',5,'FilterSize',21);

aVx = [];
aVy = [];

adrx = [];
adry = [];

for f = 1:numel(fpath3)
    fpath = fullfile(fpath1,fpath2,fpath3(f),'img');
    files = dir(fpath);
    files = files(~[files.isdir]);

    for t = 1:numel(files)
        t

        switch dataset
            case "mbwt"
                img = imread(fullfile(files(t).folder,files(t).name));
                l = rgb2gray(img);
                l = double(imgaussfilt(imsharpen(l,'Amount',3,'Radius',3),3));
                l = l./imgaussfilt(l,128);
                rescale(l,0,1);
                if (t == 1)
                    lp = l;
                end
        end

        flow = estimateFlow(opticFlow,l);
        meanVx = mean(flow.Vx,'all');
        meanVy = mean(flow.Vy,'all');
        h = xcorr_fft(l,lp);
        p = xcorrpeak(h);
        sz = size(l);
        drx = p(1) - sz(2)/2;
        dry = p(2) - sz(1)/2;

        aVx = [aVx; meanVx];
        aVy = [aVy; meanVy];
        adrx = [adrx; drx];
        adry = [adry; dry];

        lp = l;
        
    end
end


%%
plot(aVx,-adrx,'.','Color',[0.5 0.5 0])
hold on
plot(aVy,-adry,'.','Color',[0 0 0.5])
plot([-1 1],[-8 7],'r')
drawnow

set(gca,'xlim',[-1 1],'ylim',[-20 20],'FontSize',24);
