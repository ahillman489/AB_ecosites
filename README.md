Alberta Ecosite modelling

This repository contains the code for predictive ecosite mapping of Alberta. Multivariate random forest models are based on field data collected by the Alberta Biodiversity Monitoring Institute (ABMI), Ecological Information Systems (ECOSYS), Alberta Permanent Sample Plots (PSP), and individual researchers Alex Lanti-Traikovski (Canadian Forest Service), Michelle Knaggs (University of Alberta), and Ashley Hillman (Canadian Forest Service, University of Alberta).

Site locations for ABMI data within this repository use publically available coordinate data.

This repository contains three main files:

Site_Locations.csv - datasheet containing site-level coordinates, moisture regime, and nutrient regime for all sample sites.

Ecosites_climate - code for models that contain climate data (Mean Annual Temperature and Mean Annual Precipitation)

Ecosites_no_climate - code for models that do not contain climate data.

Input rasters can be found here:
https://drive.google.com/drive/folders/19L6JIxkNSqpMXXf1GasTzE5jAMkNWNeg?usp=drive_link

All scripts are written in R
