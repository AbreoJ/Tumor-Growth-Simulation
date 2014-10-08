function Visual3Dtumor3()

global N len celltype

cell=celltype;

cell(find(cell<0.95))=0;
cell(find(cell>=0.95))=1;

%========================== visualize tumor =============================
kkk=reshape(cell*3,N,N,N);

for i=1:N        %???????3D?????????????3Dtumor??????????????
      kkk(:,:,i)=kkk(:,:,i)';
end

fig=figure;

phandles = contourslice(kkk,[],[],[101],1);

daspect([1,1,1])
lightangle(45,30); 

close(fig)

figure

kkk = smooth3(kkk,'gaussian');   %Gaussian smooth filter

% kkk = smooth3(kkk);   %Gaussian smooth filter


daspect([1,1,1])
lightangle(45,30); 

hiso = patch(isosurface(kkk,1),...
'FaceColor',[1,.75,.65],...
'EdgeColor','None');

set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50)

axis([0 200 0 200 0 200])
grid on 
box on 
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Tumor and Angiogenesis')

alpha(0.4)
view(72,30)

box off

end