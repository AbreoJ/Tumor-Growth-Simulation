function Celldivide_new(s)

% This function is for tumor cell proliferation.

global celltype cell_energy stackvalue index_bias stackcount pres

stackcount=0;

nec=1e-4;

prob=celltype(s+index_bias); 
index1=(find(prob>0));
index2=(find(prob==0));
prob(index1)=0;
prob(index2)=1;

    if isempty(find(prob==1))
        prob=ones(1,length(index_bias));
        scalar0=(celltype(s+index_bias)~=nec);
        prob=prob.*scalar0;  
        
        if ~isempty(find(prob==1))
            
            prob=prob/norm(prob);
            prob=prob/sum(prob);
            prob=cumsum(prob);

            prob=[0 prob];
            mov=rand(1);

            for i=1:length(index_bias)
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

            if celltype(s+index_bias(gro))==0
               celltype(s+index_bias(gro))=1;     
               cell_energy(s)=15;
               cell_age(s+index_bias(gro))=15;
            else
               celltype(s+index_bias(gro))=1;     
               cell_energy(s)=15;
               old_cell_energy=cell_energy(s+index_bias(gro));
               cell_energy(s+index_bias(gro))=15;
               stackvalue=[s+index_bias(gro), old_cell_energy];
               cellmove()
            end
            
        end
    else
        
        scalar=pres(s+index_bias); 
        scalar=scalar-pres(s);
        index1=(find(scalar>=0));
        index2=(find(scalar<0));
        scalar(index1)=1;
        scalar(index2)=0;

        prob=prob.*scalar;
        
        if isempty(find(prob==1))

            index0=find(celltype(s+index_bias)==0);

            gro_index=randperm(length(index0));
            gro=index0(gro_index(1));

            celltype(s+index_bias(gro))=1;     
            cell_energy(s)=15;
            cell_energy(s+index_bias(gro))=15;

        else
            prob=prob/norm(prob);
            prob=prob/sum(prob);
            prob=cumsum(prob);

            prob=[0 prob];
            mov=rand(1);

            for i=1:length(prob)
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

            celltype(s+index_bias(gro))=1;     
            cell_energy(s)=15;
            cell_energy(s+index_bias(gro))=15;
        end
  
    end
end


function cellmove()    % Cells are pushed to move outwards. 

%This function is used for tumor cell movement when inner cells proliferate.
global celltype N slen cell_energy stackvalue stackcount 

stackcount=stackcount+1;

nec=1e-4;

if stackcount<=450    % To avoid stack overflow

    s=stackvalue(1);
    old_cell_energy=stackvalue(2);

    position=[N,1,slen,-N,-1,-slen];

    prob=celltype(s+position); 
    index1=(find(prob>0));
    index2=(find(prob==0));
    prob(index1)=0;
    prob(index2)=1;

    if isempty(find(prob==1))    
        prob=ones(1,length(position));
        scalar0=(celltype(s+position)~=nec);
        prob=prob.*scalar0;   
    end
     
    if ~isempty(find(prob==1))
        
        prob=prob/norm(prob);
        prob=prob/sum(prob);
        prob=cumsum(prob);

        prob=[0 prob];
        mov=rand(1);

        for i=1:length(position)
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

        if celltype(s+position(gro))==0
           celltype(s+position(gro))=1;
           cell_energy(s+position(gro))=old_cell_energy;
        else
           celltype(s+position(gro))=1;     
           old_cell_energy_move=cell_energy(s+position(gro));   
           cell_energy(s+position(gro))=old_cell_energy; 
           stackvalue=[s+position(gro), old_cell_energy_move];
           cellmove()
        end

    end
end

end



