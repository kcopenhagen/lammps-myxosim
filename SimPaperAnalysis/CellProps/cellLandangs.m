%% Matt's WT.
fpaths = '/scratch/gpfs/kc32/matt-monolayers/mb-wt/2022-11-20/';
runs = ["trial1"; "trial2"; "trial3"];
dxspl = 0.4;

eds = 0:0.02:0.5;
Leds = 0:20;
mbwtaL = [];
mbwtaangs = [];
isoangs = [];
XYcal = 0.072;
test_n = 0.01; % percent of cells per frame.
for r = 1:numel(runs)
    fpath = fullfile(fpaths,runs(r),'seg');
    files = dir(fpath);
    files = files(~[files.isdir]);
    
    for fi = 1:numel(files)
        fname = fullfile(files(fi).folder,files(fi).name);
        [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(fname);
        
        f = memmapfile(fname, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
        
        %%
        tmp = f.Data.d;
        cmsks = tmp>0.5;
        CC = bwconncomp(cmsks);
        R = regionprops(CC,'PixelIdxList','Area');
        R = R([R.Area]>42);

        Rtest = R(rand(size(R))<test_n);

        for i = 1:numel(Rtest)
            imt = zeros(size(cmsks));
            imt(Rtest(i).PixelIdxList) = 1;
            [Rspl, L] = mask2spline(imt,dxspl/XYcal);
            Rtest(i).splx = Rspl(:,1)';
            Rtest(i).sply = Rspl(:,2)';
            Rtest(i).L = L;
        end
%% Checks for isolated cells.
        % for i = 1:numel(R)
        % 
        %     Rt = R;
        %     Rt(i) = [];
        %     ixs = R(i).splx;
        %     iys = R(i).sply;
        %     oxs = [Rt.splx];
        %     oys = [Rt.sply];
        %     oxs = [oxs 1:size(imt,2) ...
        %         size(imt,2)*ones(1,size(imt,1)) ...
        %         1:size(imt,2) ...
        %         ones(1,size(imt,1))];
        %     oys = [oys size(imt,1)*ones(1,size(imt,2)) ...
        %         1:size(imt,1) ...
        %         ones(1,size(imt,2)) ...
        %         1:size(imt,1)];
        %     dx = ixs - oxs';
        %     dy = iys - oys';
        %     dr = sqrt(dx.^2 + dy.^2);
        %     if (min(dr,[],'all') > 5/XYcal)
        %         R(i).iso = true;
        %     else
        %         R(i).iso = false;
        %     end
        % end

        for i = 1:numel(Rtest)
            
            mbwtaL = [mbwtaL; Rtest(i).L];
    
            dx1 = Rtest(i).splx(1:end-2) - Rtest(i).splx(2:end-1);
            dx2 = Rtest(i).splx(2:end-1) - Rtest(i).splx(3:end);
            dy1 = Rtest(i).sply(1:end-2) - Rtest(i).sply(2:end-1);
            dy2 = Rtest(i).sply(2:end-1) - Rtest(i).sply(3:end);
            dr1 = sqrt(dx1.^2 + dy1.^2);
            dr2 = sqrt(dx2.^2 + dy2.^2);
            dot12 = dx1.*dx2 + dy1.*dy2;
    
            angs = acos(dot12./(dr1.*dr2));
            % if R(i).iso
            %     isoangs = [isoangs; abs(angs)'];
            % end
            mbwtaangs = [mbwtaangs; abs(angs)'];
        end

        % disp(sprintf('%s\n%d Isolated cells',fname,sum([R.iso])))
    end
end
% Matt's pilA
% fpaths = '/scratch/gpfs/kc32/matt-monolayers/mb-pila/2022-11-13/';
% runs = ["trial1"; "trial2"; "trial3"];
fpaths = '/scratch/gpfs/kc32/matt-monolayers/mb-pila/2022-06-23/';
runs = ["trial1"; "trial2"; "trial3"; "trial4"];

eds = 0:0.02:0.5;
Leds = 0:20;
mbpilaaL = [];
mbpilaaangs = [];
XYcal = 0.072;

for r = 1:numel(runs)
    fpath = fullfile(fpaths,runs(r),'seg');
    files = dir(fpath);
    files = files(~[files.isdir]);
    
    for fi = 1:numel(files)
        fname = fullfile(files(fi).folder,files(fi).name);
        [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(fname);
        
        f = memmapfile(fname, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
        
        %%
        tmp = f.Data.d;
        tmpt = tmp(:,:,1).*(1-tmp(:,:,2));
        cmsks = tmpt>0.1;
        CC = bwconncomp(cmsks);
        Rt = regionprops(CC,'PixelIdxList','Area');
        R = Rt([Rt.Area]>42);

        Rtest = R(rand(size(R))<test_n);
        for i = 1:numel(Rtest)
            imt = zeros(size(cmsks));
            imt(Rtest(i).PixelIdxList) = 1;
            [Rspl, L] = mask2spline(imt,dxspl/XYcal);
            mbpilaaL = [mbpilaaL; L];
    
            dx1 = Rspl(1:end-2,1) - Rspl(2:end-1,1);
            dx2 = Rspl(2:end-1,1) - Rspl(3:end,1);
            dy1 = Rspl(1:end-2,2) - Rspl(2:end-1,2);
            dy2 = Rspl(2:end-1,2) - Rspl(3:end,2);
            dr1 = sqrt(dx1.^2 + dy1.^2);
            dr2 = sqrt(dx2.^2 + dy2.^2);
            dot12 = dx1.*dx2 + dy1.*dy2;
    
            angs = acos(dot12./(dr1.*dr2));
            mbpilaaangs = [mbpilaaangs; abs(angs)];
        end
    end
end
% My WT.
XYcal = 0.072;
fpath = '/tigress/kc32/train/msk';
files = dir(fpath);
files = files(~[files.isdir]);
gd = arrayfun(@(x) x.name(1) ~= '.',files);
files = files(gd);
kcaL = [];
kcaangs = [];
kcdens = [];

for t = 1:numel(files)
    im = load(fullfile(files(t).folder,files(t).name));

    cmsks = im.msk2cls == 1;
    CC = bwconncomp(cmsks);
    R = regionprops(CC,'PixelIdxList','Area');
    R = R([R.Area]>42);
    Rtest = R(rand(size(R))<test_n);
    kcdens = [kcdens; numel(R)/(size(cmsks,1)*size(cmsks,2)*XYcal*XYcal)];
    axs = [];
    ays = [];
    for i = 1:numel(Rtest)
        imt = zeros(size(cmsks));
        imt(Rtest(i).PixelIdxList) = 1;
        [Rspl, L] = mask2spline(imt,dxspl/XYcal);
        axs = [axs; Rspl(:,1)];
        ays = [ays; Rspl(:,2)];
        kcaL = [kcaL; L];

        dx1 = Rspl(1:end-2,1) - Rspl(2:end-1,1);
        dx2 = Rspl(2:end-1,1) - Rspl(3:end,1);
        dy1 = Rspl(1:end-2,2) - Rspl(2:end-1,2);
        dy2 = Rspl(2:end-1,2) - Rspl(3:end,2);
        dr1 = sqrt(dx1.^2 + dy1.^2);
        dr2 = sqrt(dx2.^2 + dy2.^2);
        dot12 = dx1.*dx2 + dy1.*dy2;

        angs = acos(dot12./(dr1.*dr2));
        kcaangs = [kcaangs; abs(angs)];

    end
end

% Simulations

eps = 5;
d = 0.7;
fp = "PoissonRevs";
l = 7;
rho = 0.30;
v = 5;
Kagar = 300;
Kstiffs = [5 10 15 20 25 30 35 40 45 50 60 70 80 90 ...
    100 120 140 160 180 200 250 300 400 500];
rev = 8;
runs = 0;
res = 1389;
props = [];
test_n = 0.02;

for k = 1:numel(Kstiffs)
    Kstiff = Kstiffs(k);
    aLsim = [];
    aangssim = [];
    for r = 1:numel(runs)
        run = runs(r);
        fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,run,'PoissonRevs');
        [files, boxSize] = getframes(fpath);
        for t = 1:100:numel(files)
            bds = loadsimdata(fullfile(files(t).folder,files(t).name));
            aids = unique([bds.mol]);
            aids = aids(rand(size(aids))<test_n);
            
            for id = 1:numel(aids)
                ccell = bds([bds.mol]==aids(id));
                xs = [ccell.xs];
                ys = [ccell.ys];
                cspl = [xs' ys'];
                msk = spline2mask(cspl,boxSize,res,2);
                cspl2 = mask2spline(msk,dxspl/(boxSize/res),21);
                xxs = cspl2(:,1);
                yys = cspl2(:,2);

                % [~,bdinds] = sort([ccell.id]);
                % ccell = ccell(bdinds);

                % xxs = [ccell.xs];
                % yys = [ccell.ys];

                % cim = zeros(res);
                % inds = sub2ind(size(cim),ceil(yys*res/boxSize),ceil(xxs*res/boxSize));
                % cim(inds) = 1;

                dx = diff(xxs);
                dx(dx>boxSize/2) = dx(dx>boxSize/2) - boxSize;
                dx(dx<-boxSize/2) = dx(dx<-boxSize/2) + boxSize;
                dy = diff(yys);
                dy(dy>boxSize/2) = dy(dy>boxSize/2) - boxSize;
                dy(dy<-boxSize/2) = dy(dy<-boxSize/2) + boxSize;
                dr = sqrt(dx.^2 + dy.^2);
                
                aLsim = [aLsim; sum(dr)];

                dot12 = dx(1:end-1).*dx(2:end) + dy(1:end-1).*dy(2:end);

                angs = abs(acos(dot12./(dr(1:end-1).*dr(2:end))));
                aangssim = [aangssim; angs];
            end
        end
    end

    propt.Kspr = Kstiffs(k);
    propt.aLs = aLsim;
    propt.aangs = aangssim;
    props = [props; propt];
end

sims = props;

%
save(fullfile('~/SimPaperAnalysis/Figures/Data','aangaLswsim20230830.mat'),'sims','kcaL','kcaangs',...
    'mbwtaL','mbwtaangs','mbpilaaangs','mbpilaaL');

%% Plot lengths

Leds = 0:2:20;
xs = (Leds(2)-Leds(1))/2:(Leds(2)-Leds(1)):Leds(end);
mbpila = histcounts(mbpilaaL*0.076,Leds,'Normalization','pdf');
mbwt = histcounts(mbwtaL*0.076,Leds,'Normalization','pdf');
kc = histcounts(kcaL*0.133,Leds,'Normalization','pdf');
sim = histcounts(sims(8).aLs,Leds,'Normalization','pdf');
plot(xs,mbwt,'r','LineWidth',2);
hold on
plot(xs,mbpila,'g','LineWidth',2);
plot(xs,kc,'m','LineWidth',2);
plot(xs,sim,'b','LineWidth',2);

set(gca,'FontSize',24);
xlabel('Length ({\mu}m)')
ylabel('PDF')
%% Plot angles
angeds = 0:0.04:0.5;
xs = (angeds(2)-angeds(1))/2:(angeds(2)-angeds(1)):angeds(end);
mbpila = histcounts(mbpilaaangs,angeds,'Normalization','pdf');
mbwt = histcounts(mbwtaangs,angeds,'Normalization','pdf');
kc = histcounts(kcaangs,angeds,'Normalization','pdf');

simcols = sky(max([sims.Kspr]));

for s = 1:numel(sims)
    cts = histcounts(sims(s).aangs,angeds,'Normalization','pdf');
    plot(xs,cts,'Color',simcols(sims(s).Kspr,:),'LineWidth',2);
    hold on
end

%plot(xs,mbwt,'r','LineWidth',2);
%plot(xs,mbpila,'g','LineWidth',2);
%plot(xs,kc,'m','LineWidth',2);

set(gca,'FontSize',24);
xlabel('Bead angles (rad)')
ylabel('PDF');
colormap(sky)
c = colorbar;
c.YTickLabel = [min([sims.Kspr]) round(mean([min([sims.Kspr]) max([sims.Kspr])])) max([sims.Kspr])];
c.Label.String = "K_{stiff}";
%% Isolated cells vs. all cells.
% 
% eds = 0:0.02:0.5;
% xs = (eds(1) + (eds(2)-eds(1))/2):(eds(2)-eds(1)):eds(end);
% 
% % Just grouped cells.
% ysa = histcounts(mbwtaangs,eds);
% ysi = histcounts(isoangs,eds);
% ys = (ysa - ysi);
% ys = ys./(sum(ys)*(eds(2)-eds(1)));
% plot(xs,ys,'Color',[1 0.5 0],'LineWidth',2);
% hold on
% % All cells.
% ys = histcounts(mbwtaangs,eds,'Normalization','pdf');
% %plot(xs,ys,'r','LineWidth',2);
% 
% % Just isolated cells.
% ys = histcounts(isoangs,eds,'Normalization','pdf');
% plot(xs,ys,'Color',[1 0 0.5],'LineWidth',2);
% 
% % Full monolayer.
% ys = histcounts(mbpilaaangs,eds,'Normalization','pdf');
% plot(xs,ys,'g','LineWidth',2)

%%

angeds = 0:0.04:0.5;
xs = (angeds(2)-angeds(1))/2:(angeds(2)-angeds(1)):angeds(end);
mbpila = histcounts(mbpilaaangs,angeds,'Normalization','pdf');
mbwt = histcounts(mbwtaangs,angeds,'Normalization','pdf');
kc = histcounts(kcaangs,angeds,'Normalization','pdf');

kcers = [];
mbpilaers = [];
mbwters = [];

xxx = [];
for k = 1:numel(props)
    xxx = [xxx; props(k).Kspr];
    ys = histcounts(props(k).aangs,angeds,'Normalization','pdf');
    
    kcer = ((ys(kc>0) - kc(kc>0))./kc(kc>0)).^2;
    
    mbpilaer = ((ys - mbpila)./mbpila).^2;
    mbwter = ((ys - mbwt)./mbwt).^2;
    kcers = [kcers; mean(kcer,'omitnan')];
    mbpilaers = [mbpilaers; mean(mbpilaer,'omitnan')];
    mbwters = [mbwters; mean(mbwter,'omitnan')];
end

plot(xxx,kcers,'LineWidth',2,'Color','m')
hold on
plot(xxx,mbpilaers,'g','LineWidth',2);
plot(xxx,mbwters,'r','LineWidth',2);
set(gca,'FontSize',24,'YScale','log','XScale','log');
