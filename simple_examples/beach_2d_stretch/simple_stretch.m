clear all
dxmin=0.5;
dxmax=2.0;
dymax=2.0;

% grid 1, 2m constant size
m1=360;
n1=250;
dx1=2.0;
dy1=2.0;
x1=[0:m1-1]*dx1;
y1=[0:n1-1]*dy1;

[X1, Y1]=meshgrid(x1,y1);

xs=400.0;
depmax=10.0;
depmin=-2.0;
dep1=zeros(n1,m1)+depmax;
for j=1:n1
    for i=1:m1
        if x1(i)>xs
        dep1(j,i)=depmax+(x1(i)-xs)*(depmin-depmax)/(x1(m1)-xs);
        end
    end
end

% grid 2 stretched

m=500;
n=250;
dep=zeros(n,m);

m_st=100;





    for i=1:m_st
        dx(i)=dxmax;
    end
    for i=m_st+1:m
        dx(i)=dxmax+(dxmin-dxmax)*(i-m_st)/(m-m_st);
    end

x(1)=0.0;
for i=1:m-1
    x(i+1)=x(i)+dx(i);
end

y=[0:n-1]*dymax;

dy=(1:n)*0.0+dymax;

[DX,DY]=meshgrid(dx,dy);
[X, Y]=meshgrid(x,y);


dep=griddata(X1,Y1,dep1,X,Y);

xx=x';
yy=y';

save -ASCII depth.txt dep
save -ASCII depth1.txt dep1
save -ASCII dx.txt DX
save -ASCII dy.txt DY
save -ASCII x.txt xx
save -ASCII y.txt yy


    
