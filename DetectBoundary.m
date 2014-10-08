function boundary=DetectBoundary(nod3xyz,len)

% This function is designed for detecting computational domain boundary

geom.a=find(nod3xyz(:,1)==0);      %left   boundary
geom.b=find(nod3xyz(:,1)==len);    %right  boundary
geom.c=find(nod3xyz(:,2)==0);      %front  boundary
geom.d=find(nod3xyz(:,2)==len);    %back   boundary
geom.e=find(nod3xyz(:,3)==0);      %bottom boundary
geom.f=find(nod3xyz(:,3)==len);    %top    boundary

geom.c=setdiff(geom.c,[geom.a;geom.b]);
geom.d=setdiff(geom.d,[geom.a;geom.b]);
geom.e=setdiff(geom.e,[geom.a;geom.b;geom.c;geom.d]);
geom.f=setdiff(geom.f,[geom.a;geom.b;geom.c;geom.d]);

boundary=[geom.a;geom.b;geom.c;geom.d;geom.e;geom.f];

end