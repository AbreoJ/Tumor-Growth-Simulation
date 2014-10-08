function Spreadhotpoint(magnitude)

% This function is performed to generate a bunch of sprouting stimulus points spread in the 3D computational domain, the density of those sprouting
% stimulus points are calculated according to the concentration of VEGF at each given grid. 

global TAF nod3xyz N len slen wlen hotpoint

hotpoint=zeros(1,wlen);      % Reset hotpoint

% posi=magnitude*1.2.^log(TAF);

posi=magnitude*1.3.^log(TAF);

index=find(rand(wlen,1)<posi);

hotpoint(index)=1;

% scatter3(nod3xyz(index,1)*20,nod3xyz(index,2)*20,nod3xyz(index,3)*20,'filled','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','g'); 
% axis([0 200 0 200 0 200])

end
