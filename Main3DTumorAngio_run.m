
%%================ Multiscale Tumor Growth and Angiogenesis Model ===============%%
% 
% This is the main function to run a simulation. Click on the "Run" botton
% or in Command Window type in: >> Main2D_noDrug

clear all
close all
clc
tic 

%=============================    Start  ==================================
% Create a folder to save computational results
disp('==================  Simulate Tumor Growth and Anigogeneis ==================')
flnm=['TumorGrowth_Results'];           
filepath=cd;
fil=fullfile(filepath,flnm);
subfold1=strcat(filepath,'\TumorGrowth_Results\Figures');
subfold2=strcat(filepath,'\TumorGrowth_Results\Data');

if ~exist(fil)
    mkdir(fil);
end

if ~exist(subfold1)
    mkdir(subfold1);
end

if ~exist(subfold2)
    mkdir(subfold2);
end

%===================    Parameter Globalization   ===================
global len             % Virtual 3D computational domain length, unit: mm
global N               % grid number
global slen
global wlen

global nod3xyz         % wlen*3 Matrix: Matrix"Tumor Cell Locaiton in 3D Spaces"

global nutr            % 3D Matrix: Nutrient
global pres            % 3D Matrix: Pressure
global TAF             % 3D Matrix: TAF
global waste           
global activity        % Cell activity
global celltype        % 3D Matrix: Cell type
global cell_energy     % 3D Matrix: "Tumor Cell energy for proliferation"
global drug

global vess            % 3D Cell matrix: "Endothelial Cell Information"
global vess_tag        % 3D Matrix denotes "Endothelial Cell"   1:normal EC,  0.95: Tip EC
global vess_age        % 3D Matrix: "Vessel Cell Age" 

global hotpoint        % Vessel branching hotpoint

global index_bias      % 3D Matrix: offsets toward neighboring grids from one certain grid point
global stackvalue      % Vector: transfer parameter values for stack 
global stackcount      % Stack operation counter 
global vessgrowth_flag   
global branchrecord    % Record the vessel cells that have branched, we assure the vessel cells that have branched won't branch again.  
global sprout_index   
global nec


%===================    Create grid   ===================
% The domain of simulations is a cube sized 1 cm^3
% We assume the size of each tumor cell is 50 micrometers
% Hence, in each cube length there are 200 grids (or max 200 cells) 
len     = 10;           
N       = 201;             
slen    = N*N;           
wlen    = N*N*N;         


%===================    Create Vectors   ===================
%%% Create vectors of initial values of all variables
nutr        = ones(1,wlen);     % Normalized nutrient distribution, original nutrient density is set to be 1. 
waste       = ones(1,wlen);     
TAF         = zeros(1,wlen);
drug        = zeros(1,wlen);


%===================    Initialization other variables  ===================
pres        = zeros(1,wlen);     
activity    = zeros(1,wlen);  % Tumor cell activity should be zero at initial stage
celltype    = zeros(1,wlen);
cell_energy = zeros(1,wlen);   

vess        = cell(1,wlen);          
vess_tag    = zeros(1,wlen);
vess_age    = zeros(1,wlen);

hotpoint        = zeros(1,wlen);
branchrecord    = zeros(1,wlen);

stackvalue      = [0 0];
vessgrowth_flag = 0;
stackcount      = 0;

xn          = linspace(0,len,N);
[Y X Z]     = meshgrid(xn);  % Note: Y X Z should be in this order
nod3xyz     = [X(:),Y(:),Z(:)];

nec     = 1e-4;  % Indicator of necrosis cells

clear X Y Z fil flnm xn filepath subfold1 subfold2

boundary = DetectBoundary(nod3xyz,len);  % Boundary identification
InitCancerCell();                        % Initialization cancer cell number & position
hold on
Vasculature3D();                         % Draw artifical vasculature
hold off     


%=======================   Initial vessel cells  ===========================
% Initial sprouting points located along artificial vasculature.
X = [100 135 172 183 187 181 157 123 86 52 24 19 45 18 10];
Y = [186 179 163 123 87 47 25 15 14 23 44 148 180 77 115];
Z = repmat(101,1,15);

varalist = {'TAF', 'activenumber', 'act', 'activity', 'cell_energy', ...
    'celltype', 'necrosisnumber', 'quiescentnumber', 'pres', 'nutr', ...
    'v_age', 'vess', 'vess_age', 'vess_tag', 'waste', 'drug'};

vindex = zeros(size(X));

for i = 1:length(X)
    xx = round(X(i)/20/0.05);
    yy = round(Y(i)/20/0.05);
    zz = round(Z(i)/20/0.05);
    
    vindex(i) = xx+yy*N+zz*slen;    
end

hold on

scatter3(nod3xyz(vindex,1)*20,nod3xyz(vindex,2)*20,nod3xyz(vindex,3)*20,'g*');

for j = 1:length(vindex)
    s = vindex(j);
    vess{s}.count   = 0;
    vess{s}.pare    = [];
    vess{s}.son     = [];
    vess{s}.direct  = [];
end

vess_tag(vindex) = 0.95;
vess_age(vindex) = 1;

vindexsave = vindex;
old_vindexsave = vindex;

clear vindex i j k X Y Z xx yy zz

%======================     Update System Status    =======================
totDays     = 60;  % Total simulation days (two months)
tau         = totDays*24*3600;  % Unit: second
L           = 1e-2;  % Unit: m, computational domain length: 10 mm
k           = 1/(len*(N-1));  % Normalized time step = 0.0005.     2000 iterations
h           = 1/(N-1);  % Spatial step size for numerical method - central difference is used for convection term

c0      = 4.3e-4;  % Standard TAF concentration, unit: kg/m^3,  || Refer to "significant expression of vascular endothelial growth factor/vascular permeability factor in mouse ascites tumors" 
n0      = 8.4;  % Standard nutrient concentration, unit: mol/m^3,  || Refer to "HEMOGLOBIN-BASED OXYGEN CARRIER HEMODILUTION AND BRAIN OXYGENATION"
w0      = 10.5;  % Standard waste concentration
d0      = 2.13;  % Standard drug concentration, unit: mol/m3


% Nutrient parameters
Dn          = 8e-14;  % Nutrient diffusion rate, unit: m^2/s [Wang & Li, 1998]
rho_n0      = 6.8e-4;  % Vessel nutrient supply rate, unit: mol/(m^3*sec) [Wang & Li, 1998]
lambda_n0   = 3.0e-5;  % Nutrient consumption rate [Vaupel et al, 1987]

% Waste parameters
Dw         = 4e-14;  % Carbon dioxide diffusion coefficient, unit: m^2/s [estimated] 
rho_w0     = 1e-5;  % Carbon dioxide secretion rate, unit: mol/(m^3s) [estimated]     
lambda_w0  = 2.5e-5;  % Carbon dioxide consumption rate, unit: ml/(cm^3s) [estimated]

% TAF parameters
Dc          = 1.2e-13;  % TAF diffusion coefficient [estimated] 
rho_c0      = 2e-9;  % TAF secretion rate [estimated]
lambda_c0   = 0;  % VEGF natural decay rate, assumed to be very small [estimated]
 
% Drug parameters
Dd              = 1.5e-14;  % Drug diffusion coefficient [estimated] 
lambda_d0       = 2.5e-7;  % Drug consumption rate [estimated]
lambda_decay    = 1.0e-8;  % Drug decay rate [estimated]

% Pressure
pres_scale  = 1;  % Tumor intersitial presssure: 1~60 mmHg, by adjusting this parameter, we can therefore investigate the growth patterns of different types of tumor, low pressure, high pressure.
                  % Our hypothesis is that the interstitial pressure inside solid tumor will give rise to different morphologies. Bigger interstitial pressure gives rise to dendritic tumor.
cap_pres    = 30;  % unit: mmHg
p0          = 60*pres_scale;  % unit: mmHg

Kp          = 4.5e-15;  % Hydraulic conductivity of the interstitium, unit: cm^2/mmHg-sec
u0          = 5e-6;
pv          = cap_pres/p0;  % Capillary/vascular pressure, unit: (20) mmHg 
sigma0      = 0.15;
amplitude0  = 0.08*pres_scale;  % Gaussian function amplitude

Lp          = 2.8e-9;  % Hydraulic conductivity of the microvascular wall, unit: m/mmHg-sec [Baxter & Jain, 1989]
sigmaT      = 0.82;  % Average osmotic reflection coefficient in tumor [Baxter & Jain, 1989]  
sigmaD      = 0.1;
Pvp         = 1.49e-9;  % Vascular permeability coefficient, unit:  m/sec [estimated]
piV         = 0.3546/p0;  % Osmotic pressure of the drug, unit: mmHg
piI         = 0.2667/p0;  % Osmotic pressure of the interstitial fluid, unit: mmHg
 
k_AR2       = 500;


%%%%%%%%%%  index_bias is used to evaluate tumor cell density  %%%%%%%%%%
kn=1;
index_bias0=[];
index_bias=[];
seq=-kn:kn;
for i=-kn:kn
   index_bias0=[index_bias0,seq+N*i];
end
for j=-kn:kn
   index_bias=[index_bias,index_bias0+slen*j];
end

clear i j

%%%%%%%%%%  pres_bias is used to calculate tumor pressure  %%%%%
knn=20;
pres_bias0=[];
pres_bias=[];
seq=-knn:knn;
for i=-knn:knn
   pres_bias0=[pres_bias0,seq+N*i];
end
for j=-knn:knn
   pres_bias=[pres_bias,pres_bias0+slen*j];  
end

clear index_bias0 pres_bias0 totDays i j kn knn seq 

t1=setdiff(1:wlen,boundary);      % Find inner points by removing boundary cells

figure

sprout_index=[];


Cd=1.0*[zeros(1,1000),ones(1,1000)];    % Blood drug concentration  

% Uncomment these two lines if you want to start the simulation from a
% certain day, e.q., day 40th. But for this, you need to have an input of 
% simulation at day 40.
% filename=['./TumorGrowth_Results/Data'];
% load([filename,'/DataDay' num2str(40) '.mat'])

% 1 day = 33 iterations
% 40 days = 1321 iterations
% 60 days = 1981 iterations

calit = 1981;

for iteration=1:calit  % (Drug treatment: Day40->Day60) (Day0-->Day40 without Drug)
    
    if iteration >= 1321
        ddf = 1.0;
    else
        ddf = 0.0;
    end
    
    iteration  

    activenumber(iteration)=length(find(celltype==1));         % Number of living cell
    quiescentnumber(iteration)=length(find(celltype==0.95));   % Number of quiescent cell
    necrosisnumber(iteration)=length(find(celltype==nec));     % Number of necrotic cell
    
    disp('Cell number of different phenotypes:')
    fprintf('\n %s : %d       %s : %d        %s : %d\n\n', 'Living cell:',activenumber(iteration),'Quiescent cell:', quiescentnumber(iteration), 'Necrosis cell:', necrosisnumber(iteration));
    plot(iteration,sum([activenumber(iteration)  quiescentnumber(iteration)  necrosisnumber(iteration)]),'b*'), drawnow;
    hold on 
    plot(iteration,activenumber(iteration),'ro'),    drawnow;
    plot(iteration,quiescentnumber(iteration),'gs'), drawnow;
    plot(iteration,necrosisnumber(iteration),'k^'),  drawnow;
    
    t=reshape(TAF,N,N,N);
    n=reshape(nutr,N,N,N);
    v=reshape(vess_tag,N,N,N);
    c=reshape(celltype,N,N,N);
    w=reshape(waste,N,N,N);
    v_age=reshape(vess_age,N,N,N);  
    act=reshape(activity,N,N,N);
    d=reshape(drug,N,N,N);
    
    v_rad = v_age./(k_AR2+v_age);
    
    %%%%%%%%%%  calculate pressure  %%%%%%%%%%
    disp('Calculating pressure...') 
    
   
    %=== Calculate every 20th step ===%
    if mod(iteration-1,20)==0       
        pres=zeros(1,wlen);        % Pressure field needs to be recalculated each iteration.    
        cindex=find(celltype~=0);  % Pls note: We are not interested in cells, since it is not the cells cause intersitial pressure to its surrounding tissue. We are just interested in each grind point occupied by tumor cells which indicate the tumor region.  
        vindex=find(vess_tag>0);   % All vessel cells 
        
        % Calculate CTP
        for tt=1:length(cindex) 
           s=cindex(tt);
           
           % Density as defined by Equation (3c)
           density  = sum(celltype(s+index_bias)>0)/length(index_bias);
           
           % Euclidean distance in Equation (2), (4a)
           dist=sqrt(sum((nod3xyz(s+pres_bias,:)-repmat(nod3xyz(s,:),length(pres_bias),1)).^2,2));
           
           % Alpha in Equation (4b) is assumed to be constant
           amplitude=amplitude0;
           
           % Sigma as defined by Equation (4c)
           sigma=sigma0*density^2/(density^2.+0.5^2)+0.05;
           
           % Pressure based on Gaussian function as in Equation (4a)
           pres(s+pres_bias)=pres(s+pres_bias) + amplitude*(exp(-dist'.^2./(2*sigma^2)));
        end
        
        % Calculate VTP
        for ss = 1:length(vindex)
            vind = vindex(ss);
            
            % Tumor cell density around vessel cell, here we assume vessel
            % only causes additional pressure inside tumor due to membrane
            % effect [Ref]
            density  = sum(celltype(vind+index_bias)>0)/length(index_bias);
            
            % Equation(4b), alpha0 set to be 0.01
            amplitude = 0.01*density^2/(density^2.+0.5^2); 

            %Euclidean distance in Equation (4a)
            dist=sqrt(sum((nod3xyz(vind+pres_bias,:)-repmat(nod3xyz(vind,:),length(pres_bias),1)).^2,2));
            
            % Equation (4c)
            sigma=sigma0*density^2/(density^2.+0.5^2)+0.05;
            
            % Add vessel cell induced pressure 
            pres(vind+pres_bias)=pres(vind+pres_bias) + amplitude*(exp(-dist'.^2./(2*sigma^2)));
        end
        
        %%%%%%%%%%  Compute Pressure Gradient  %%%%%%%%%%
        gradp=zeros(3,wlen); 
 
        pres=pres/p0;  
        
        gradp(1,t1) = (pres(t1+1) - pres(t1-1))/ (2*h); 
        gradp(2,t1) = (pres(t1+N) - pres(t1-N))/ (2*h);
        gradp(3,t1) = (pres(t1+slen)-pres(t1-slen)) /(2*h);
        
        ux=-Kp*p0/(u0*L)*reshape(gradp(1,:),N,N,N);
        uy=-Kp*p0/(u0*L)*reshape(gradp(2,:),N,N,N);
        uz=-Kp*p0/(u0*L)*reshape(gradp(3,:),N,N,N);

        ux=reshape(ux,N,N,N);
        uy=reshape(uy,N,N,N);
        uz=reshape(uz,N,N,N);       
    end
    

   
    p = reshape(pres,N,N,N);
    weight = (cap_pres-p.*p0)./cap_pres;    %60-->30mmHg, 30-->25mmHg
    weight(find(weight<0)) = 0;
    

    
    %%%%%%%%%%  Nutrient Equation  %%%%%%%%%%
    disp('Calculating nutrient...')
        
    c1=1-6*Dn*tau*k/(L^2*h^2); 
    c2=Dn*tau*k/(L^2*h^2)-tau*k*u0/L*ux(2:N-1,2:N-1,2:N-1)/(2*h); 
    c3=Dn*tau*k/(L^2*h^2)-tau*k*u0/L*uy(2:N-1,2:N-1,2:N-1)/(2*h); 
    c4=Dn*tau*k/(L^2*h^2)-tau*k*u0/L*uz(2:N-1,2:N-1,2:N-1)/(2*h); 
    c5=Dn*tau*k/(L^2*h^2)+tau*k*u0/L*ux(2:N-1,2:N-1,2:N-1)/(2*h); 
    c6=Dn*tau*k/(L^2*h^2)+tau*k*u0/L*uy(2:N-1,2:N-1,2:N-1)/(2*h); 
    c7=Dn*tau*k/(L^2*h^2)+tau*k*u0/L*uz(2:N-1,2:N-1,2:N-1)/(2*h); 

    rho_n = rho_n0*v_rad.*weight; 
    lambda_n = lambda_n0*act(2:N-1,2:N-1,2:N-1);
    
    rho_nf = k*tau/n0*rho_n(2:N-1,2:N-1,2:N-1).*v(2:N-1,2:N-1,2:N-1);
    lambda_nf = k*tau/n0*lambda_n.*c(2:N-1,2:N-1,2:N-1);
    
    n(2:N-1,2:N-1,2:N-1)= c1.*n(2:N-1,2:N-1,2:N-1) + c2.*n(3:N,2:N-1,2:N-1) + c3.*n(2:N-1,3:N,2:N-1) ...
        + c4.*n(2:N-1,2:N-1,3:N) + c5.*n(1:N-2,2:N-1,2:N-1) + c6.*n(2:N-1,1:N-2,2:N-1) ...
        + c7.*n(2:N-1,2:N-1,1:N-2) ...
        + rho_nf - lambda_nf;   
                         
    temp_n=n;
    temp_n(boundary)=1;
    nutr=temp_n(:)';
    
    
    %%%%%%%%%%  Metabolic Waste Equation  %%%%%%%%%%
    disp('Calculating Metabolic Waste...')
    
    a1=1-6*Dw*tau*k/(L^2*h^2); 
    a2=Dw*tau*k/(L^2*h^2)-tau*k*u0/L*ux(2:N-1,2:N-1,2:N-1)/(2*h); 
    a3=Dw*tau*k/(L^2*h^2)-tau*k*u0/L*uy(2:N-1,2:N-1,2:N-1)/(2*h); 
    a4=Dw*tau*k/(L^2*h^2)-tau*k*u0/L*uz(2:N-1,2:N-1,2:N-1)/(2*h); 
    a5=Dw*tau*k/(L^2*h^2)+tau*k*u0/L*ux(2:N-1,2:N-1,2:N-1)/(2*h); 
    a6=Dw*tau*k/(L^2*h^2)+tau*k*u0/L*uy(2:N-1,2:N-1,2:N-1)/(2*h); 
    a7=Dw*tau*k/(L^2*h^2)+tau*k*u0/L*uz(2:N-1,2:N-1,2:N-1)/(2*h); 

    w_new=w;
    w_new(find(w_new<1))=1;
    act=n./(n+1)*2.*exp(-(w_new-1).^4/0.2);
    
    rho_w = rho_w0*act(2:N-1,2:N-1,2:N-1);
    lambda_w = lambda_w0*weight;
    
    rho_wf = k*tau/w0*rho_w.*c(2:N-1,2:N-1,2:N-1);
    lambda_wf = k*tau/w0*lambda_w(2:N-1,2:N-1,2:N-1).*v(2:N-1,2:N-1,2:N-1).*v_rad(2:N-1,2:N-1,2:N-1);
    
    w(2:N-1,2:N-1,2:N-1) = a1.*w(2:N-1,2:N-1,2:N-1) + a2.*w(3:N,2:N-1,2:N-1) + a3.*w(2:N-1,3:N,2:N-1) ...
        + a4.*w(2:N-1,2:N-1,3:N) + a5.*w(1:N-2,2:N-1,2:N-1) + a6.*w(2:N-1,1:N-2,2:N-1) ...
        + a7.*w(2:N-1,2:N-1,1:N-2) ...
        + rho_wf - lambda_wf;   
                         
    temp_w=w;
    temp_w(boundary)=1;
    waste=temp_w(:)';
      
    
    %%%%%%%%%%  TAF Equation  %%%%%%%%%%
    disp('Calculating VEGF...')
          
    s1=1-6*Dc*tau*k/(L^2*h^2);
    s2=Dc*tau*k/(L^2*h^2)-tau*k*u0/L*ux(2:N-1,2:N-1,2:N-1)/(2*h); 
    s3=Dc*tau*k/(L^2*h^2)-tau*k*u0/L*uy(2:N-1,2:N-1,2:N-1)/(2*h); 
    s4=Dc*tau*k/(L^2*h^2)-tau*k*u0/L*uz(2:N-1,2:N-1,2:N-1)/(2*h); 
    s5=Dc*tau*k/(L^2*h^2)+tau*k*u0/L*ux(2:N-1,2:N-1,2:N-1)/(2*h); 
    s6=Dc*tau*k/(L^2*h^2)+tau*k*u0/L*uy(2:N-1,2:N-1,2:N-1)/(2*h); 
    s7=Dc*tau*k/(L^2*h^2)+tau*k*u0/L*uz(2:N-1,2:N-1,2:N-1)/(2*h); 
    
    n_new=n(2:N-1,2:N-1,2:N-1);  n_new(find(n_new>1))=1;
    
    rho_c = rho_c0*(1-n_new);
    lambda_c = lambda_c0*v_rad;
    
    rho_cf = k*tau/c0*rho_c.*c(2:N-1,2:N-1,2:N-1);
    lambda_cf = k*tau/c0*lambda_c(2:N-1,2:N-1,2:N-1).*v(2:N-1,2:N-1,2:N-1);
    
    t(2:N-1,2:N-1,2:N-1) = s1.*t(2:N-1,2:N-1,2:N-1) + s2.*t(3:N,2:N-1,2:N-1) + s3.*t(2:N-1,3:N,2:N-1) ...
        + s4.*t(2:N-1,2:N-1,3:N) + s5.*t(1:N-2,2:N-1,2:N-1) + s6.*t(2:N-1,1:N-2,2:N-1) ...
        + s7.*t(2:N-1,2:N-1,1:N-2) + rho_cf - lambda_cf;
                    
    t(boundary)=0;
    TAF=t(:);
    activity=act(:)';  
    
    
    %%%%%%%%%%  Drug Equation  %%%%%%%%%%
    disp('Calculating Drug...')
    
 
    v_age(find(v_age==0))=1e-8;  % To avoid NaN during calculation
    SV=250*v_rad;  % unit: /m, mature vessel supplies more nutrient compared to young vessels
    SV=SV(2:N-1,2:N-1,2:N-1);  
    
    p=reshape(pres,N,N,N);
   
    Delta_P = (pv-p(2:N-1,2:N-1,2:N-1)-sigmaT*(piV-piI));
    Delta_P(find(Delta_P<=0)) = 0;  % We hypothesize that drug won't go back into blood vessel again. 
    
    Fv = Lp*p0*SV.*Delta_P;
    Pev = Lp*p0*Delta_P*(1-sigmaD)./Pvp;
 
    r1 = 1-6*Dd*tau*k/(L^2*h^2);
    r2 = Dd*tau*k/(L^2*h^2)-tau*k*u0/L*ux(2:N-1,2:N-1,2:N-1)/(2*h); 
    r3 = Dd*tau*k/(L^2*h^2)-tau*k*u0/L*uy(2:N-1,2:N-1,2:N-1)/(2*h); 
    r4 = Dd*tau*k/(L^2*h^2)-tau*k*u0/L*uz(2:N-1,2:N-1,2:N-1)/(2*h); 
    r5 = Dd*tau*k/(L^2*h^2)+tau*k*u0/L*ux(2:N-1,2:N-1,2:N-1)/(2*h); 
    r6 = Dd*tau*k/(L^2*h^2)+tau*k*u0/L*uy(2:N-1,2:N-1,2:N-1)/(2*h); 
    r7 = Dd*tau*k/(L^2*h^2)+tau*k*u0/L*uz(2:N-1,2:N-1,2:N-1)/(2*h); 
    
    c_new=c(2:N-1,2:N-1,2:N-1);  c_new(find(c_new>1))=0;  c_new(find(c_new==0))=1;
    
    rho_d = Fv.*(1-sigmaD)*d0.*Cd(iteration) + Pvp*SV*d0.*(Cd(iteration)-d(2:N-1,2:N-1,2:N-1)).* Pev./((exp(Pev)-1)+Inf);
    lambda_d1 = lambda_d0*act(2:N-1,2:N-1,2:N-1);
    lambda_d2 = lambda_decay*d(2:N-1,2:N-1,2:N-1);
    
    rho_df = k*tau/d0*rho_d.*v(2:N-1,2:N-1,2:N-1);
    lambda_d1f = k*tau/d0*lambda_d1.*c(2:N-1,2:N-1,2:N-1).*d(2:N-1,2:N-1,2:N-1).*d(2:N-1,2:N-1,2:N-1);
    lambda_d2f = k*tau/d0*lambda_d2.*c_new;
    
    d(2:N-1,2:N-1,2:N-1)= r1.*d(2:N-1,2:N-1,2:N-1) + r2.*d(3:N,2:N-1,2:N-1) + r3.*d(2:N-1,3:N,2:N-1) ...
        + r4.*d(2:N-1,2:N-1,3:N) + r5.*d(1:N-2,2:N-1,2:N-1) + r6.*d(2:N-1,1:N-2,2:N-1) ...
        + r7.*d(2:N-1,2:N-1,1:N-2) ...
        + rho_df ...
        - lambda_d1f ...
        - lambda_d2f;
       
                   
    d(boundary)=0;
    drug=ddf*d(:)';
    
%     clear c c1 c2 c3 c4 c5 c6 c7 t n p v TAF_weight temp_n Lp SV sigma sigmaT Pvp pv pi ci cv Pev n_new rho_n
    clear s1 s2 s3 s4 s5 s6 s7 r1 r2 r3 r4 r5 r6 r7 a1 a2 a3 a4 a5 a6 a7
    

    %  ====================  Tumor Growth model  =========================
    disp('Calculating tumor growth...')  
  
    cellindex=find(celltype>nec);                          % Living <--> Quiescent, Quiescent --> Necrotic
    test_act=activity(cellindex);
    test_cell_energy=cell_energy(cellindex);
    celltype(cellindex(find(test_act>=0.5)))=1;            % Active tumor cells
    celltype(cellindex(find(test_act<0.5)))=0.95;          % Quiescent tumor cells
    celltype(cellindex(find(test_cell_energy<=0)))=nec;    % Necrosis Cells
    
    % Cell proliferation
    div_index=find(celltype==0.95);
    div_index=div_index(randperm(length(div_index)));
    
    kill=10;
    
    for ss=1:length(div_index)
        s=div_index(ss);      
        cell_energy(s)=cell_energy(s)-0.1; 
    end

    div_index=find(celltype==1);
    div_index=div_index(randperm(length(div_index)));
    
    for ss=1:length(div_index)
        s=div_index(ss);
        prolif_energy=30; 
              
        if cell_energy(s)>=prolif_energy     
           Celldivide_new(s);
        else
           cell_energy(s) = cell_energy(s) + activity(s) - activity(s)./(activity(s)+1) ...
               - drug(s)*activity(s)./(activity(s)+1)*kill;  
        end
    end
    
    
    % 3D Angiogenesis 
    startime=400;
    
    if iteration>=startime         
          disp('Calculating 3D angiogenesis ...')
           
          vessgrowth_flag=1;
          start=startime;

          if iteration==startime
              tip_index=vindexsave;
          end

         vess_age(find(vess_tag>0))=vess_age(find(vess_tag>0))+1;   % vess_age++ in each iteration
         tip_index=union(find(vess_tag==0.95),sprout_index);        % Including the new sprouting points lying on preexisting blood vessel that were touched by floating hotpoints and other known tip cells.                                            
         tip_index=tip_index(randperm(length(tip_index)));          % Randperm tip_index

         for i=1:length(tip_index)
             s=tip_index(i); 
             Angiogenesis3D(s);
         end

        Spreadhotpoint(0.3e-3);  % To spread vessel branching hotpoints in 3D computational domain
        sprout_index=Sproutcheck();  % Check if the vessel cell would sprout.(Hot points touching blood vessel casuses branching)
        
    end
     
   %%%%%%%%%%  Save figures & data and draw tumor  %%%%%%%%%%
   
   if mod(iteration-1,33)==0
        num_day=round((iteration-1)/33); 
        Visual3Dtumor()
        filename=cd;
        saveas(gca,[filename,'\TumorGrowth_Results\Figures\TumorDay' num2str(num_day) '_' num2str(1) '.fig'])
        close
        Visual3Dtumor3()
        filename=cd;
        saveas(gca,[filename,'\TumorGrowth_Results\Figures\TumorDay' num2str(num_day) '_' num2str(2) '.fig'])
        close
   end
   
   hold off
   
   if mod(iteration-1,33)==0
        num_day=round((iteration-1)/33); 
        filename=cd;
        save([filename,'\TumorGrowth_Results\Data','\DataDay' num2str(num_day) '.mat'],varalist{:})
   end
   
   hold on
    
end

hold off

totaltime=toc

Tumorvessel3D();  % draw 3D tumor and vessel 

Visual3Dtumor()
hold on
box off
alpha(0.6)

xlabel('x')   
ylabel('y')  
zlabel('z')
title('nutrient, z-crosssection')

hold on
Draw3Dvessel()




