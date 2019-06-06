clear all
dxmin=0.5;
dxmax=2.0;
dxsponge=100.0;
dymax=2.0;

m=550;
n=250;

m_sponge_st=200;
m_sponge=100;
m_st=300;

for i=1:m_sponge
        dx(i)=dxsponge;
end

for i=m_sponge+1:m_sponge_st
        dx(i)=dxsponge+(dxmax-dxsponge)*(i-m_sponge)/(m_sponge_st-m_sponge);
end
for i=m_sponge_st+1:m_st
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


xs=900.0;
depmax=10.0;
depmin=-2.0;
dep=zeros(n,m)+depmax;

%for j=1:n
 %   for i=1:m
  %      if x(i)>xs
   %     dep(j,i)=depmax+(x(i)-xs)*(depmin-depmax)/(x(m)-xs);
    %    end
   % end
% end


xx=x';
yy=y';

save -ASCII depth.txt dep
save -ASCII dx.txt DX
save -ASCII dy.txt DY
save -ASCII x.txt xx
save -ASCII y.txt yy


    
