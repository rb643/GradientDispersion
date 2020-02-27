
## Cam CAN analyses of gradient dispersion
This repo contains all code used in *_"Dispersion of functional gradients across the lifespan"_* by Bethlehem & Paquola et al.


### Main code and folders
**./PostProcessing/camcan_functional_gradients.m**    
This function runs the main gradient analyses including:
1. Gradient construction    
2. Linear effects within individual gradient    
3. Projection of gradients into 3D gradient space   
4. Quantification of within and between network dispersion    
5. Generate surface projections for figures     
6. Saves all output to csv for follow analyses and figure plotting

**./PostProcessing/camcan_mediation.R**
Runs the mediation analyses  

**./PostProcessing/camcan_moderation.R**
Runs the moderation analyses

**./PostProcessing/connectomise_holdout.m**
Generates connectomes from surface files

**./Figures_code/...**
Contains code used to generate paper figures

**./PreProcessing/...**
Contains all code used to pre-process data

**./Utilities/...**
Contains various matlab functions, colourmaps and utiltities for plotting. Also includes a setup script to setup folder structure locally (just change your relevant folder paths :)

**./Data/...**
All data used for analyses including connectomes and output generated during analyses, excluding the Conte69 surface maps (due to space constraints)

#### System
Tested on linux 18.04/16.0, using Matlab v17-19

#### Dependencies
- Matlab (for gradient construction and analyses)
- SurfStat (for surface based projection, analyses and visualisation)
- MICA Tools (e.g. Matlab tools and Surfstat add-ons): https://github.com/MICA-MNI/micasoft    
- R (for figures and mediation and moderation analyses)
