clear all
close all
fdir1='../beach_2d_stretch/output/';
fdir2='../beach_2d_constant/output/';
fmodel='../beach_2d_stretch/';

eta=load([fdir2 'eta_00001']);

x1=load([fmodel 'x.txt']);
y1=load([fmodel 'y.txt']);

[n,m]=size(eta);
dx2=2.0;
dy2=2.0;
x2=[0:m-1]*dx2;
y2=[0:n-1]*dy2;


nfile=[10 30 50 70];

wid=8;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf


for num=1:length(nfile)
    
fnum=sprintf('%.5d',nfile(num));
eta1=load([fdir1 'eta_' fnum]);
mask1=load([fdir1 'mask_' fnum]);

eta1(mask1==0)=NaN;


figure(num)
subplot(1,2, 1)

pcolor(x1,y1,eta1),shading flat
hold on
caxis([-0.5 1.5])
title(['Stretched, t = ' num2str(num*2) ' sec '])
axis([0 700 0 500])

ylabel(' y (m) ')
xlabel(' x (m) ')

subplot(1,2, 2)

eta2=load([fdir2 'eta_' fnum]);
mask2=load([fdir2 'mask_' fnum]);
eta2(mask2==0)=NaN;

pcolor(x2,y2,eta2),shading flat
hold on
caxis([-0.5 1.5])
title(['Constant, t = ' num2str(num*2) ' sec '])
axis([0 700 0 500])
ylabel(' y (m) ')
xlabel(' x (m) ')

name=num2str(num*2);

print('-djpeg100', ['eta_',name, 'min.jpg']);

end
%print -djpeg eta_inlet_shoal_irr.jpg