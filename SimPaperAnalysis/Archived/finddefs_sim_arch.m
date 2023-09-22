function adefs = finddefs_sim(bds,boxSize,varargin)
%% Input - bds, boxSize. 
% Optional input arguments - XYcal (grid size), 
% Sr (Order boxsize), 
% Dr (director smoothing), 
% Scutoff (maximum S for defects), 
% loopr (loop for charge calc.

    XYcal = 0.2;
    Sr = 1;
    Dr = 1;
    Scutoff = 0.99;
    loopr = 0.8;

    if mod(nargin,2)~=0
        disp("Invalid input format");
    end
    while numel(varargin)>0
        switch varargin{1}
            case 'XYcal'
                XYcal = varargin{2};
                varargin(1:2) = [];
            case 'Sr'
                Sr = varargin{2};
                varargin(1:2) = [];
            case 'Dr'
                Dr = varargin{2};
                varargin(1:2) = [];
            case 'Scutoff'
                Scutoff = varargin{2};
                varargin(1:2) = [];
            case 'loopr'
                loopr = varargin{2};
                varargin(1:2) = [];
            otherwise
                'Invalid input'
        
        end
        
    end
    
    dirf = dfield_sim(bds,boxSize,XYcal,Dr);
    S = nemorderfield_sim(dirf,Sr);
    S(S>Scutoff) = 1;

    S = imgaussfilt(S,2);
    Smins = imregionalmin(S,8);
    %posdefs = Smins.*(1-S)>(1-Scutoff);
    
    [defy, defx] = find(Smins);

    [re,rimsz] = makeringelement(round(loopr/XYcal));

    padsz = round(round(rimsz/2)+1);

    dirfpad = padarray(dirf,[padsz padsz],'circular');
    defxpad = defx + padsz;
    defypad = defy + padsz;
    qs = zeros(size(defx));
    d = [];
    for i = 1:numel(defx)
        dirt = dirfpad(defypad(i)-(rimsz-1)/2:defypad(i)+(rimsz-1)/2,...
            defxpad(i)-(rimsz-1)/2:defxpad(i)+(rimsz-1)/2);
        angs = dirt(re);
        dangs = diff(angs);
        dangs(dangs>pi/2) = dangs(dangs>pi/2) - pi;
        dangs(dangs<-pi/2) = dangs(dangs<-pi/2) + pi;
        qs(i) = sum(dangs)/(2*pi);
        
        if (abs(qs(i)-0.5)<0.1)
            nx = cos(dirt);
            nxnx = nx.*nx;
            ny = sin(dirt);
            nyny = ny.*ny;
            nxny = nx.*ny;
            [dnxnxdx, ~] = gradient(nxnx);
            [dnxnydx, dnxnydy] = gradient(nxny);
            [~, dnynydy] = gradient(nyny);
            dnxnxdx = dnxnxdx((rimsz+1)/2,(rimsz+1)/2);
            dnxnydx = dnxnydx((rimsz+1)/2,(rimsz+1)/2);
            dnxnydy = dnxnydy((rimsz+1)/2,(rimsz+1)/2);
            dnynydy = dnynydy((rimsz+1)/2,(rimsz+1)/2);
            dx = dnxnxdx+dnxnydy;
            dy = dnxnydx+dnynydy;
            dr = sqrt(dx^2+dy^2);
            d = [d; [dx/dr dy/dr]];
        elseif abs(qs(i)+0.5) < 0.1
            
            nx = cos(-dirt);
            nxnx = nx.*nx;
            ny = sin(-dirt);
            nyny = ny.*ny;
            nxny = nx.*ny;
            
            [dnxnxdx, ~] = gradient(nxnx);
            [dnxnydx, dnxnydy] = gradient(nxny);
            [~, dnynydy] = gradient(nyny);
            dnxnxdx = dnxnxdx((rimsz+1)/2,(rimsz+1)/2);
            dnxnydx = dnxnydx((rimsz+1)/2,(rimsz+1)/2);
            dnxnydy = dnxnydy((rimsz+1)/2,(rimsz+1)/2);
            dnynydy = dnynydy((rimsz+1)/2,(rimsz+1)/2);
            
            psiprime = atan2(dnxnydx+dnynydy,dnxnxdx+dnxnydy);
            d = [d; [cos(-psiprime/3) sin(-psiprime/3)]];
        else
            d = [d; NaN NaN];
        end
    end

    del = zeros(size(defx));
    del(abs(qs)<0.1) = 1;

    defx(del==1) = [];
    defy(del==1) = [];
    qs(del==1) = [];
    d(del==1,:) = [];

    adefs = struct('x',num2cell(defx*XYcal),'y',num2cell(defy*XYcal),...
        'q',num2cell(qs),'dx',num2cell(d(:,1)),'dy',num2cell(d(:,2)));


    
end