clear all
fdir='../input_files/output/';

eta=load([fdir 'eta_00001']);

[n,m]=size(eta);
dx=2.0;
dy=2.0;
x=[0:m-1]*dx;
y=[0:n-1]*dy;

x_sponge=[0 100 100 0 0];
y_sponge=[0 0 y(end) y(end) 0];
x_wavemaker=[150 160 160 150 150];
y_wavemaker=[0 0 y(end) y(end) 0];


nfile=[10 50];
min={'20' '100'};

wid=8;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf


for num=1:length(nfile)
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);
mask=load([fdir 'mask_' fnum]);

eta(mask==0)=NaN;

subplot(1,length(nfile), num)

pcolor(x,y,eta),shading flat
hold on
caxis([-0.5 2])
title([' Time = ' min{num} ' sec '])

hold on
plot(x_sponge,y_sponge,'g--','LineWidth',2)
h1=text(50,500,'Sponge','Color','w');
set(h1, 'rotation', 90)

plot(x_wavemaker,y_wavemaker,'r-','LineWidth',2)
h2=text(180,700,'Wavemaker','Color','w');
set(h2, 'rotation', 90)

if num==1
ylabel(' y (m) ')
end

xlabel(' x (m) ')
%cbar=colorbar;
%set(get(cbar,'ylabel'),'String','\eta (m) ')

set(gcf,'Renderer','zbuffer')

end
%print -djpeg eta_inlet_shoal_irr.jpg