function  Vasculature3D() 

% Generate a 3D vasculature image based on hand draw 2D
% vasculature image named vesselimage2D.bmp. 

filepath=cd;
filepath=strcat(filepath,'\vesselimage2D.bmp');
data=imread(filepath);
image01=data(:,:,2);
image01=im2double(image01);
[N,N]=size(image01);
test=image01(:);
index=find(test<1);

% Generate a 3D curve
[y,x] = meshgrid([1:1000]);

z1 = 200*exp(-x.^2/(2*380^2));  % surface1
z2 = 300*exp(-((x-450).^2+(y-480).^2)/(2*300^2))-60;  % surface2
z3 = 200*exp(-(x-1050).^2/(2*100^2));  % surface3
z4 = 200*exp(-(y-1080).^2/(2*100^2));  % surface4
z5 = -150*exp(-((x-1000).^2+(y-1000).^2)/(2*100^2));  % surface5
z= (z1+z2+z3+z4+z5+50)/2;

% Draw 3D images after Processing
yy0=(index-mod(index,N))/N;
xx0=index-yy0*N;
zz0=test(index)+z(index);

yy=((index-mod(index,N))/N)/5;
xx=(index-yy0*N)/5;
zz=(test(index)+z(index)-42)/147*200;

index=find(xx==0);
xx(index)=[];
yy(index)=[];
zz(index)=[];

scatter3(xx,yy,zz,40,'r.')

axis([0 200 0 200 0 200])
title('3D Artificial Vasculature')
view(72,30)

end
