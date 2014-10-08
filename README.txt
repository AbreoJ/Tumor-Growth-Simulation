Matlab code for running simulations in the manuscript titled:
"Multiscale Tumor Growth and Angiogenesis Model"

The package contains the following Matlab files:
1. Main3DTumorAngio_run.m, the main file to run a simulation
2. InitCancerCell.m
3. DetectBoundary.m
4. Angiogenesis3D.m
5. Draw3Dvessel.m
6. Vasculature3D.m
7. Spreadhotpoint.m
8. Sproutcheck.m
9. Tumorvessel3D.m
10. Visual3Dtumor.m
11. Visual3Dtumor3.m
12. Celldivide_new.m

an an image file named "vesselimage2D.bmp" for vessel initialisation.

To run a simulation, open file "Main3DTumorAngio_run.m" in Matlab editor and click the "Run" botton or press F5 or from Matlab's command window type in "Main3DTumorAngio_run" (in the same directory).

All parameters in Table 1 are presented in file Main3DTumorAngio_run.m. Modify accordingly.

The simulations will create a folder "TumorGrowth_Results" and subfolders "Data" and "Figures" where .mat data and figures will be saved, respectively. To undo this, comment out lines 576-594.