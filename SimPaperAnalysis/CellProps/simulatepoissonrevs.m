N = 100;
T = 30;
l = 7.43;
dt = 0.1;
d = 0.27;
figure
p = (1/l)*dt;

ttrs = [];

for n = 1:N
    ttrs = [ttrs; struct('revs',[0])];
end

for t = 0:dt:(T-dt)
    for i = 1:numel(ttrs)
        if (rand() < p) && (abs(ttrs(i).revs(end) - t) > d)
            ttrs(i).revs = [ttrs(i).revs; t];
        end
    end
end

acrevs = [];
arevs = [];
for i = 1:numel(ttrs)
    crevs = [ttrs(i).revs; T];
    lifetimes = diff(crevs);
    if (numel(lifetimes) == 1)
        arevs = [arevs; lifetimes(1) 1 1];
    else
        arevs = [arevs; lifetimes(1) 1 0];
        arevs = [arevs; lifetimes(end) 0 1];
    end
    if (numel(lifetimes)>2)
        for t = 2:numel(lifetimes)-1
            arevs = [arevs; lifetimes(t) 0 0];
        end
    end

    acrevs = [acrevs; lifetimes];
end

eds = 0:0.5:60;
xs = (eds(2)-eds(1))/2:(eds(2)-eds(1)):eds(end);
ys = histcounts(acrevs,eds);

histogram(acrevs,eds)
hold on
fitl = fit(xs',ys','exp1');
plot(fitl)
fitl
set(gca,'FontSize',24,'xlim',[0 30])
