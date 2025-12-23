# Using cold hardiness dynamics to predict budbreak and low temperature risk in grapevines

Code and data used for predictions of cold hardiness and phenology using the NYUS.1 model, and phenology predictions using PhenoFlex model for comparison.

## NYUS.1 Analysis

The folder contains:

### Data files

- VitisPhenology_data_budbreak.csv: budbreak data with Cultivar, Year, bb.date (date of budbreak), bb.doy (budbreak day of year), Location, BBCH.scale used in observation, Observer (if available), Dataset (season in which predictions are being made, spanning two years)
- VitisPhenology_data_locations.csv: file used to add Latitude information to each location for calculation of hourly temperatures
- VitisPhenology_data_weather.csv: compiled weather data, with daily Tmin and Tmax
- VitisPhenology_data_Figure4 RMSE chilling variation.csv: data used for figure 4, for which manual changes in chilling calculation were made.

### Analysis files

- VitisPhenology_functions.r: functions that were created for cold hardiness and phenology predictions
- VitisPhenology_NYUS.1 Model_data process and figure prep.Rmd: Analysis and figure preparation code.

### How to run the analysis

1. All NYUS.1 analysis files should be in the same folder
2. Open .Rmd file
3. All data and functions are loaded, and predictions are produced. 
4. Steps within analysis are annotaded within code
5. Figures are created individually

## Supplementary PhenoFlex Analysis

This analysis is independent from that of NYUS.1.

The folder contains:

### Data files

- VitisPhenology_PhenoFlex_data_budbreak.csv
- VitisPhenology_PhenoFlex_data_locations.csv
- VitisPhenology_PhenoFlex_data_weather.csv

### Analysis file

- VitisPhenology_PhenoFlex_data process and figure prep.Rmd

### How to run the analysis

1. All PhenoFlex analysis files should be in the same folder
2. Open .Rmd file
3. All data files are loaded, and predictions are produced. 
4. Figures are created individually
