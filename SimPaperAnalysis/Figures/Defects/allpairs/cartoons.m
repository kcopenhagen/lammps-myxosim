%% np cartoon

plot(0,0,'b.','MarkerSize',30);
hold on
plot([0 1],[0 0],'b','LineWidth',3);
plot([0 cos(2*pi/3)],[0 sin(2*pi/3)],'b','LineWidth',3)
plot([0 cos(4*pi/3)],[0 sin(4*pi/3)],'b','LineWidth',3);
axis equal
axis off

thetas = 0:((2*pi)/12):2*pi-0.01;
phis = -thetas*2+pi;

phis(phis>pi) = phis(phis>pi)-2*pi;
phis(phis>pi) = phis(phis>pi)-2*pi;
phis(phis<-pi) = phis(phis<-pi)+2*pi;
phis(phis<-pi) = phis(phis<-pi)+2*pi;
cmapinds = round(255*(phis+pi)/(2*pi))+1;

colorcet('C9')
cmap = get(gca,'Colormap');

for thi = 1:numel(thetas)
    plot(3*cos(thetas(thi)),3*sin(thetas(thi)),'.','MarkerSize',30,...
        'Color',cmap(cmapinds(thi),:))
    plot([3*cos(thetas(thi)) 3*cos(thetas(thi)) + cos(phis(thi))],...
        [3*sin(thetas(thi)) 3*sin(thetas(thi)) + sin(phis(thi))],...
        'Color',cmap(cmapinds(thi),:),'LineWidth',3);
end

set(gca,'YDir','reverse');

%% pn cartoon

plot(0,0,'r.','MarkerSize',30);
hold on
plot([0 1],[0 0],'r','LineWidth',3);
axis equal
axis off

thetas = 0:((2*pi)/12):2*pi-0.01;
phis = thetas*2+pi;

phis(phis>pi) = phis(phis>pi)-2*pi;
phis(phis>pi) = phis(phis>pi)-2*pi;
phis(phis<-pi) = phis(phis<-pi)+2*pi;
phis(phis<-pi) = phis(phis<-pi)+2*pi;
cmapinds = round(255*(phis+pi)/(2*pi))+1;

colorcet('C9')
cmap = get(gca,'Colormap');

for thi = 1:numel(thetas)
    plot(3*cos(thetas(thi)),3*sin(thetas(thi)),'.','MarkerSize',30,...
        'Color',cmap(cmapinds(thi),:))
    plot([3*cos(thetas(thi)) 3*cos(thetas(thi)) + cos(phis(thi)/3)],...
        [3*sin(thetas(thi)) 3*sin(thetas(thi)) + sin(phis(thi)/3)],...
        'Color',cmap(cmapinds(thi),:),'LineWidth',3);
    plot([3*cos(thetas(thi)) 3*cos(thetas(thi)) + cos(phis(thi)/3+2*pi/3)],...
        [3*sin(thetas(thi)) 3*sin(thetas(thi)) + sin(phis(thi)/3+2*pi/3)],...
        'Color',cmap(cmapinds(thi),:),'LineWidth',3);
    plot([3*cos(thetas(thi)) 3*cos(thetas(thi)) + cos(phis(thi)/3+4*pi/3)],...
        [3*sin(thetas(thi)) 3*sin(thetas(thi)) + sin(phis(thi)/3+4*pi/3)],...
        'Color',cmap(cmapinds(thi),:),'LineWidth',3);
end

set(gca,'YDir','reverse');