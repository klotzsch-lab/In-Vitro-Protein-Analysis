# In Vitro Protein Analysis

The scripts in this repository are used to analyze in vitro protein data.

## Droplet Analysis

To calculate droplet characteristics like area, intensity, feret diameter etc., the python based ImageJ script "In_vitro_droplet_analysis.py" is used. This will result in a .csv file with all characteristics.

Afterwards "In_vitro_droplet_analysis.R" can be used to create graphs and for statistics in RStudio.

## FRAP Analysis

"FRAPfit.py" is used to analyze FRAP experiment data of in vitro droplets. Before using "FRAPfit.py" .csv files with intensities of FRAP trajectories are needed.