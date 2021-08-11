
README.md for scripts to accompany the **"Exploring trade-offs in mixed fisheries by integrating fleet dynamics into multispecies size-spectrum models"** paper by Novaglio et al. 

The master script RunModel.Rmd is structured into three main steps: 

**Step 1** loads libraries and data, and calls all helper codes needed to run model and scenarios. 

**Step 2** runs the calibrated coupled size-spectrum and fleet dynamics model representing Australia’s South East Shark and Scalefish Fishery and produces Fig 3. 

**Step 3** sets up scenarios, runs projections and calculates ecological and socio-economic indicators detailed in the paper. It produces Figs 4 and 5.  

** Data is saved in the /data folder, which is part of this repository 

*** Figures are not saved withing the repository for space limits

*** project.r includes project functions for the size based modelling package mizer here modified to integrate fleet dynamics 
