function Angiogenesis3D(s)

% This function is used for generating 3D neovessels in 3D computational
% domain.

global vess pres vess_tag vess_age hotpoint index_bias branchrecord sprout_index

% if vess{s}.count>=(8+20*pres(s))   % Here would not necessarily be 4, which might change according to location conditions, like pressure and tissue density.    
if vess{s}.count>=(10*((1/10)^(1/30))^(30-pres(s)*60)) 
    
   gro=direction(s,index_bias);   
   
   if ~isempty(gro) 
       if vess_tag(s)==1   

               vess{s}.son=s+index_bias(gro);      % vess{s}.son=[] when it was created, now it has one descendant, I add the son location. 
               vess_tag(s+index_bias(gro))=0.95; 
%                vess_age(s+index_bias(gro))=vess_age(s);      % New endothelial cell :  cellage=0; cellage++ in each iteration.
               if vess_age(s) > 1
                    vess_age(s+index_bias(gro))=1;      % New endothelial cell :  cellage=0; cellage++ in each iteration.
               else
                    vess_age(s+index_bias(gro)) = vess_age(s)/2;
               end
               vess{s+index_bias(gro)}.count=0;    % New endothelial cell :  count=0; count++ in each iteration.
               vess{s+index_bias(gro)}.pare=s;     % Add new EC parent location.
               vess{s+index_bias(gro)}.son=[];     % Add new EC son location, [] here.
               vess{s+index_bias(gro)}.direct=index_bias(gro);    % Record movement direction, the descendant cells will mimic their parents' behavior, used in "direction.m"
               
               branchrecord(s)=1;

               sprout_index=setdiff(sprout_index,s);

        else if vess_tag(s)==0.95

                if hotpoint(s)==0       

                   vess{s}.son=s+index_bias(gro);      % vess{s}.son=[] when it was created, now it has one descendant, I add the son location. 
                   vess_tag(s)=1;
                   vess_tag(s+index_bias(gro))=0.95;  
                   if vess_age(s) > 1
                        vess_age(s+index_bias(gro))=1;      % New endothelial cell :  cellage=0; cellage++ in each iteration.
                   else 
                        vess_age(s+index_bias(gro)) = vess_age(s)/2;
                   end
                   vess{s+index_bias(gro)}.count=0;    % New endothelial cell :  count=0; count++ in each iteration.
                   vess{s+index_bias(gro)}.pare=s;     % Add new EC parent location.
                   vess{s+index_bias(gro)}.son=[];     % Add new EC son location, [] here.
                   vess{s+index_bias(gro)}.direct=index_bias(gro);
                   
                else    % Touching the branching hotpoint, EC starts to branch.  

                   vess{s}.son=s+index_bias(gro);      % vess{s}.son=[] when it was created, now it has one descendant, I add the son location. 
                   vess_tag(s)=1;
                   vess_tag(s+index_bias(gro))=0.95;   
                   if vess_age(s) > 1
                         vess_age(s+index_bias(gro))=1;      % New endothelial cell :  cellage=0; cellage++ in each iteration.
                   else 
                         vess_age(s+index_bias(gro)) = vess_age(s)/2;
                   end
                   vess{s+index_bias(gro)}.count=0;    % New endothelial cell :  count=0; count++ in each iteration.
                   vess{s+index_bias(gro)}.pare=s;     % Add new EC parent location.
                   vess{s+index_bias(gro)}.son=[];     % Add new EC son location, [] here. 
                   vess{s+index_bias(gro)}.direct=index_bias(gro);
                   
                   index_bias2=setdiff(index_bias,index_bias(gro));
                   
                   gro=direction(s,index_bias2);
                   
                   if ~isempty(gro)

                       vess{s}.son=s+index_bias2(gro);      % vess{s}.son=[] when it was created, now it has one descendant, I add the son location. 
                       vess_tag(s)=1;
                       vess_tag(s+index_bias2(gro))=0.95;   
                       if vess_age(s)>1
                            vess_age(s+index_bias2(gro))=1;      % New endothelial cell :  cellage=0; cellage++ in each iteration.
                       else
                            vess_age(s+index_bias2(gro)) = vess_age(s)/2;
                       end
                       vess{s+index_bias2(gro)}.count=0;    % New endothelial cell :  count=0; count++ in each iteration.
                       vess{s+index_bias2(gro)}.pare=s;     % Add new EC parent location.
                       vess{s+index_bias2(gro)}.son=[];     % Add new EC son location, [] here.
                       vess{s+index_bias2(gro)}.direct=index_bias2(gro);  
                   
                   end
                   
                   branchrecord(s)=1;

                end
            end
       end
   end
   
   vess{s}.count=vess{s}.count+1;
   
else
    
    vess{s}.count=vess{s}.count+1;
    
end


end




function gro=direction(s,index_bias2)

% This function is used to generate one direction for dividing. 

global N slen TAF pres vess vess_tag celltype 

    nec=1e-4; 
    
    direct_pare=vess{s}.direct;
       
    if isempty(vess{s}.direct)
        scalar=ones(length(index_bias2),1);
    else
        [A B C]=intersect((direct_pare+[1 -1 N -N slen -slen]),index_bias2);
        scalar=zeros(length(index_bias2),1);
        scalar(C)=1;
    end
    
    prob0=TAF(s+index_bias2)-repmat(TAF(s),length(index_bias2),1);
    prob=prob0;
    index1=(find(prob>=0));      %%%%%%%% !!!!!!!!!! Before it is '>'
    index2=(find(prob<0));
    prob(index1)=1;
    prob(index2)=0;
    
    prob=prob0.*prob;     
    
    scalar2=(celltype(s+index_bias2)~=nec)';  % The tumor blood vessel won't grow into necrosis core

    prob=prob.*(vess_tag(s+index_bias2)==0)'.*scalar.*scalar2;
    
    gro=[];

    if ~isempty(find(prob>0))
        
        if (rand(1)<pres(s))   % More likely to enter this calculation when approaching the tumor 
            
            %%%% Pressure Calculation
            pres0=pres(s+index_bias2)'-repmat(pres(s),length(index_bias2),1);
            prob=pres0.*(prob==1);
            prob=repmat(norm(prob),length(index_bias2),1).*(prob~=0)-prob;
            
        end
              
        prob=prob/norm(prob);
        prob=prob/sum(prob);
        prob=cumsum(prob);

        prob=[0; prob];
        mov=rand(1);

        for i=1:length(index_bias2)
            if mov==0
               gro=1;
               break
            else
                if mov>prob(i) && mov<=prob(i+1)
                   gro=i; 
                   break
                end
            end
        end
    else
        gro=[];
    end
    
    
    
end
    

