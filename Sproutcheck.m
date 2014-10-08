function sprout_index=Sproutcheck()

% This function is performed to check whether the pre-existing ECs touch
% branching hotpoints. 

global vess_tag hotpoint nod3xyz branchrecord

sprout_index=[];

vess_index=find(vess_tag==1);

vess_index=setdiff(vess_index, find(branchrecord==1)); 

% index=find(hotpoint==1);
% num=length(index);
% 
% for i=1:length(vess_index)
%     s=vess_index(i);
%     minimum=min(sum((repmat(nod3xyz(s,:),num,1)-nod3xyz(index,:)).^2));
%     if minimum<residual
%         sprout_index=[sprout_index s];
%     end
% end

sprout_index=vess_index(find(hotpoint(vess_index)==1));


end

