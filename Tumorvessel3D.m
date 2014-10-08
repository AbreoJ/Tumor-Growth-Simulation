function Tumorvessel3D()

% This function is designed to visualize 3D solid tumor and blood vessels. 

Visual3Dtumor();

hold on 

Draw3Dvessel();

view(75.0, 16.0)

box off

end