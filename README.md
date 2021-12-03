
README.md for scripts to accompany the **"Exploring trade-offs in mixed fisheries by integrating fleet dynamics into multispecies size-spectrum models"** paper by Novaglio et al. 2021, Journal of Applied Ecology

For more information, model use and interpretation, please contact Camilla Novaglio at  camilla.novaglio@gmail.com

The master script RunModel.Rmd is structured into three main steps: 

**Step 1** loads libraries and data, and calls all helper codes needed to run model and scenarios. 

**Step 2** runs the calibrated baseline model and coupled size-spectrum and fleet dynamics model representing Australia’s South East Shark and Scalefish Fishery and produces Fig 2 and Fig 3. 

**Step 3** sets up scenarios, runs projections and calculates ecological and socio-economic indicators detailed in the paper. It produces Figs 4 and 5, as well as additional figures in Supplementary Information.  

** Data is saved in the /data folder, which is part of this repository 

*** Figures are not saved within the repository for space limits

*** project.r includes project functions from the size based modelling package mizer here modified to integrate fleet dynamics

 
