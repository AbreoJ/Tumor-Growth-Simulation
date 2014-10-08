function InitCancerCell()

% This function is designed to generate initial cancer cells located in the
% center of computational domain. 

%%================ 3D Tumor Growth and Angiogenesis Model ===============%%
%                                                                         %
%            Zhou's Lab, Methodist Hospital Research Institute            %      
%            Weill Cornell Medical College, Cornell University            %
%            Terry L.Tang            Copyright(C) Feb.13, 2012            %
%=========================================================================%


global celltype cell_energy N len slen nod3xyz activity 

h=len/(N-1);                    % Length of spatial step

num_init_cell=5;              % Initial tumor cell number
init_size=0.25;               % Initial tumor size, radius, considering len=10, usually init_size=0.5, 1.0 etc. 


index0=100*slen+100*N+101;      % Initial cancer cells location
init_cell=nod3xyz(index0,:);    % Search for center initial tumor cell

celltype(index0)=1;
cell_energy(index0)=15;         %  Half of the total energy (30) required for cell dividing

Mat_radius = [0.1793    0.0928    0.1199    0.1402];
Mat_thetal1 = [0.9848    0.2331    0.5748    0.8741];
Mat_thetal2 = [0.8167    0.9745    0.5811    0.2300];

if num_init_cell>1
    for i=1:num_init_cell-1
%         radius=rand(1)*init_size;
%         theta1=rand(1); theta2=rand(1);
        radius=Mat_radius(i);
        theta1=Mat_thetal1(i);
        theta2=Mat_thetal2(i);
        z=round(radius*sin(2*pi*theta1)/h);
        x=round(radius*cos(2*pi*theta1)*cos(2*pi*theta2)/h);
        y=round(radius*cos(2*pi*theta1)*sin(2*pi*theta2)/h);
        index=index0+x+y*N+z*slen;
        init_cell=[init_cell;nod3xyz(index,:)];
        celltype(index)=1;
        cell_energy(index)=15;
        activity(index)=1;
    end
end

init_cell=init_cell*20;

figure(1)
scatter3(init_cell(:,1),init_cell(:,2),init_cell(:,3),'filled','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g'); 
title('3Dimensional Tumor growth and Angiogenesis Model Initialization')
axis([0 200 0 200 0 200])
xlabel('x')   
ylabel('y')  
zlabel('z')
grid on

end



