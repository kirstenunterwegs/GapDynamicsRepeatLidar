# (Forest) Gap dynamic analysis with repeat lidar

This repository holds code for running the analysis of forest gap dynamics (gap formation and gap closure) basing on gap raster layers generated from repeated lidar data. For a detailed description of the purpose, methodology and results, see the following paper: 

Krüger, K., Senf, C., Jucker, T., Pflugmacher, D., Seidl, R. (2024). Gap expansion is the dominant driver of canopy openings in a temperate mountain forest landscape. Journal of Ecology. (link)

The data to reproduce the analysis can be obtained under the following link: https://zenodo.org/records/10966262

All data to reproduce the results are available, all other layers can be generated with the code provided in this repository. The Canopy Height Models underlying the analysis and respective processing scripts for the lidar data, are available from the corresponding author upon reasonable request. 


## platforms

software used for development and implementation : R version 4.3.3

## Directory structure and files:

- `data/processed/`: Any data altered from its raw format. 
- `data/processed/gaps_final`: Contains the gap layers, which were derived from Canopy Height Models (which were derived from lidar point clouds). Layers are masked to the closed forest area of the National Park core zone. Artifacts were masked out manually.
- `data/processed/gap_features`: Environmental features per individal gap in list format.
- `data/processed/gap_change`: Raster layers decoding the area of gap formtion and gap closure per time step.
- `data/processed/environment_features`: Raster layers of elevation bands, aspect and forest types within the study area.
- `data/processed/creation`: All layers generated in the analysis regarding gap formation and the classification into new and expanding gaps.
- `data/processed/closure`: All layers generated in the analysis regarding gap closure and the classification into lateral and vertical closure. Two sub-folders storing the delineation of gap boundaries.
- `data/processed/sensitivity_analysis`: Modified gap layers for a subset of the study area, to check on robustness of gap change metrics. Layers were modified one time in the minimum gap size and one time in the canopy height threshold.
- `data/results/`: Figures displaying analysis results
