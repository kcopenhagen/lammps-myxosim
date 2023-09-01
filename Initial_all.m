function Initial_all()

% N = 1843; % Maximum number of cells
rho = @n;
r_min = 1.0; %Minimum distance between beads of separate cells.
lmu = @L; % Mean length of cells
lsig = lmu/2; % Standard deviation of cell lengths
lmin = 2;
lmax = 12;
cell_w = @d;
r_eq = cell_w*0.6; %Equilibrium bead spacing (set in lammps input file).

afsig = 0.1; % Spread in active force.

boxSize_xy = 100;
boxSize_z = 50;

rng('shuffle')

xs = [0; 0];
ys = [0; 0];
zs = [-10; -10];
aids = [0; 0]; % Atom ids.
cids = []; % Cell ids.
fails = 0;
ct = 1;
totcells = 0;

while totcells < (boxSize_xy*boxSize_xy*rho) && fails < 100
    dr = 0;
    while min(dr,[],'all') < r_min
        thetat = 2*pi*rand();
        xst = lmax + (boxSize_xy - 2*lmax)*rand();
        yst = lmax + (boxSize_xy - 2*lmax)*rand();
        zst = boxSize_z/4*rand()+0.1*boxSize_z;

        lent = normrnd(lmu,lsig);
        while (lent < lmin) || (lent > lmax)
            lent = normrnd(lmu,lsig);
        end

        clen = 0;
        i = 0;
	aidt = aids(end)+1;
        xtest = [];
        ytest = [];
        ztest = [];
	aidst = [];
        while clen < lent
            xtest = [xtest; xst+r_eq*i*cos(thetat)];
            ytest = [ytest; yst+r_eq*i*sin(thetat)];
            ztest = [ztest; zst];
	    aidst = [aidst; aidt];
	    aidt = aidt+1;
            i = i+1;
            clen = clen + r_eq;
        end
        mind = boxSize_xy;
        dx = xs - xtest';
        dy = ys - ytest';
        dz = zs - ztest';
        dr = sqrt(dx.*dx+dy.*dy+dz.*dz);
        fails = fails+1;
    end
    fails = 0;
    totcells = totcells + 1;
    xs = [xs; xtest];
    ys = [ys; ytest];
    zs = [zs; ztest];
    aids = [aids; aidst];
    cids = [cids; ct*ones(size(xtest))];
    ct = ct+1;
end
xs(1:2) = [];
ys(1:2) = [];
zs(1:2) = [];
aids(1:2) = [];

ucids = unique(cids);

N = numel(xs);
Nbonds = N - numel(ucids);
Nangles = N - 2 * numel(ucids);

%% Create lammps input data file.

fID = fopen('spfr_init.txt','w');
fprintf(fID,'# Data input file for modeling myxo. \n\n');

fprintf(fID,'\t%d\tatoms\n',N);
fprintf(fID,'\t%d\tbonds\n',Nbonds);
fprintf(fID,'\t%d\tangles\n\n',Nangles);

fprintf(fID,'\t1\tatom types\n');
fprintf(fID,'\t1\tbond types\n');
fprintf(fID,'\t1\tangle types\n\n');

fprintf(fID,'\t%f\t%f\txlo xhi\n', 0, boxSize_xy);
fprintf(fID,'\t%f\t%f\tylo yhi\n\n', 0, boxSize_xy);
fprintf(fID,'\t%f\t%f\tzlo zhi\n\n', -boxSize_z, boxSize_z);

fprintf(fID,'Atoms\n\n');
% mid = aids(floor(lens(1)/2));
mid = round(mean(aids(cids==cids(1))));

af = sign(rand()*2-1)*normrnd(1.0,afsig);
i = 1;
mux = -1;
muy = 0;
muz = 0;

% Cell atoms
fprintf(fID,'\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',aids(i),mid,1,xs(i),ys(i),zs(i),mux,muy,muz,af);
for i = 2:N
    if cids(i) ~= cids(i-1)
        af = sign(rand()*2-1)*normrnd(1.0,afsig);
	mid = round(mean(aids(cids==cids(i))));
        % mid = aids(i - 1 + floor(lens(i)/2));
    end
    fprintf(fID,'\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',aids(i),mid,1,xs(i),ys(i),zs(i),mux,muy,muz,af);
end

% Cell bonds
fprintf(fID,'\nBonds\n\n');
bid = 1;
for i = 2:N
    if cids(i) == cids(i-1)
        fprintf(fID,'\t%d\t1\t%d\t%d\n',bid,i-1,i);
        bid = bid + 1;
        
    end
end

% Cell angles
fprintf(fID,'\nAngles\n\n');
aid = 1;
for i = 3:N
    if cids(i) == cids(i-2)
        fprintf(fID,'\t%d\t1\t%d\t%d\t%d\n',aid,i-2,i-1,i);
        aid = aid + 1;
    end
end

fclose(fID);

disp("Initialization file generated");
end
