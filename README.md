# (Forest) Gap dynamic analysis with repeat lidar

This repository holds code for running the analysis of forest gap dynamics (gap formation and gap closure) basing on gap layers extracted from repeated lidar data. For a detailed description of the purpose, methodology and results, see the following paper: 

Krüger, K., Senf, C., Jucker, T., Pflugmacher, D., Seidl, R. (2024). Gap expansion is the dominant driver of canopy openings in a temperate mountain forest landscape. Journal of Ecology. (link)

The gap layer generated and used for the analysis in this reserach can be obtained under the following link: (zenodo link)

## platforms

software used for development and implementation : R version 4.3.3

## Directory structure and files:

- `data/`: folder containing all input and output data as well as intermediate data layer
- `data/raw/`: Inputs and project files used in the analysis. Gap layer are available under the above stated zenodo repository.
- `data/processed/`: Any data altered from its raw format. Is dubsivided into a folder containing the gap layers (gaps_final), environmental variables (environment_features), data around gap formation (creation) and gap closure (closure).
- `data/results/`: Figures displaying analysis results
