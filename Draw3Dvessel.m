function Draw3Dvessel()

% This function is designed to illustrate 3D tumor vasculature.

global N slen vess vess_tag vess_age pres

vess_index=find(vess_tag>0);

[junk,sortindex]=sort(-vess_age(vess_index));

vess_index=vess_index(sortindex);

% vasculature3D(); 

hold on

for i=1:length(vess_index)
    
    s=vess_index(i);
    
    z0=(s-mod(s,slen))/slen;
    y0=(mod(s,slen)-mod(mod(s,slen),N))/N;
    x0=s-z0*slen-y0*N;
    
    s2=vess{s}.pare;
    
    z1=(s2-mod(s2,slen))/slen;
    y1=(mod(s2,slen)-mod(mod(s2,slen),N))/N;
    x1=s2-z1*slen-y1*N;
    
%     plot3([x0 x1], [y0 y1], [z0 z1],'-r','LineWidth',3*vess_age(s)/(vess_age(s)+250));  % I changed it from 250 --> 500   4/7/2012
    plot3([x0 x1], [y0 y1], [z0 z1], '-r','LineWidth', 3*vess_age(s)/(vess_age(s)+500));
    
%     drawnow;
    
end
    
hold off

end
    
   
    
