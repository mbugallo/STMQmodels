
# SMALL AREA ESTIMATION USING SPATIO-TEMPORAL \newline M-QUANTILE MODELS

GUIDE FOR USERS - APPLICATION TO REAL DATA

### María Bugallo¹, Domingo Morales¹, Nicola Salvati², Francesco Schirripa-Spagnolo²

### ¹Center of Operations Research, Miguel Hernández University of Elche, Spain. \newline Email of the corresponding author: mbugallo@umh.es

### ²Department of Economics and Management, University of Pisa, Italy.

## Objective

The main objective of this guide is to illustrate the application of the spatio-temporal M-quantile (STMQ) models for small area estimation (SAE) using real U.S. data. We aim to provide users with a clear workflow for data preparation, model fitting and interpretation of results, highlighting best practices and recommendations. While the datasets used are public and freely available, the focus is on demonstrating the methodology rather than distributing the data itself, ensuring that users can reproduce the analyses with the most up-to-date sources.

## Description of the real data 

This section provides a detailed overview of the data sources used in our case study. The data are organized by geographic regions (states and counties) and temporal units (years), facilitating the application of spatio-temporal modeling techniques. Details of the data sources are as follows (copied from Appendix 1 of the Supplementary Material):

1. Annual summary data from air quality measurements collected at outdoor monitoring stations across the U.S., provided by the U.S. Environmental Protection Agency (EPA), U.S. Department of the Interior. Available at: https://www.epa.gov/outdoor-air-quality-data \newline
Save the data as follows: `Target_variables/annual_aqi_by_county_YYYY.csv`, where `YYYY` represents the year (for example, 2023).

2. Small Area Income and Poverty Estimates (SAIPE) from the U.S. Census Bureau, U.S. Department of Commerce. Available at: https://www.census.gov/programs-surveys/saipe.html \newline
Save the data as follows: `Aux1_SAIPE/estYYall.csv`, where `YY` represents the year (for example, 23).

3. Regional Economic Accounts data from the U.S. Bureau of Economic Analysis, U.S. Department of Commerce. Available at: https://www.census.gov/data/datasets.html \newline
Save the data as follows: `Aux3_POP_State/ACSST1YYYY.S0101-Data.csv` (state level data) and `Aux3_POP/ACSST1YYYY.S0101-Data.csv` (county level data), where `YYYY` represents the year.

4. American Community Survey (ACS) data provided by the U.S. Census Bureau, U.S. Department of Commerce. Available at: https://www.bea.gov/data/economic-accounts/regional \newline
Save the data as follows: `Aux2/CAGDP1__ALL_AREAS_2001_2023.csv`

5. Coastal counties information from the U.S. Geological Survey (USGS), U.S. Department of the Interior. Available at: https://www.usgs.gov/the-national-map-data-delivery/ \newline
Save the data as follows: `Aux4_CoastlineCounties/coastline_counties.csv`

6. Rural-urban continuum codes from the Economic Research Service. Available at: https://www.ers.usda.gov/data-products/rural-urban-continuum-codes \newline
Save the data as follows: `Aux4_CoastlineCounties/Ruralurbancontinuumcodes2023.csv`

Data can be provided to the authors upon request; however, we recommend downloading the most up-to-date versions directly from the original sources. We prefer not to distribute the data ourselves. All datasets are public and free of charge. This study covers the period from 2016 to 2023; therefore, the years considered are 2016, 2017, 2018, 2019, 2020, 2021, 2022 and 2023.

## Installation and loading of required packages for the AQI analysis

Run the installation lines only once. After installing, keep it commented.

> install.packages(c("dplyr", "tidyr", "reshape2",  "sf", "tigris", "ggspatial",  "ggplot2", "scales", "viridis", "RColorBrewer", "MASS", "forecast", "agricolae"))

- library(dplyr)        # Data manipulation
- library(tidyr)        # Reshaping data
- library(reshape2)     # Melting and casting data frames

- library(sf)           # Spatial data manipulation
- library(tigris)       # U.S. shapefiles
- library(ggspatial)    # Map scales and north arrows

- library(ggplot2)      # Plotting
- library(scales)       # Axis scales
- library(viridis)      # Color scales
- library(RColorBrewer) # Color palettes

- library(MASS)         # Statistical functions
- library(forecast)     # Time series
- library(agricolae)    # Experimental designs

For easier use of the code, it is recommended to use the RStudio graphical user interface.

## R Code structure and description

This section provides an overview of the structure and organization of our R code. The goal is to give readers a clear understanding of the workflow, from data preparation to model fitting, prediction and visualization.   The project was developed using R version 4.3.1 and RStudio version 2024.04.1+748. RStudio requires R version 3.6.0 or higher. To the best of our knowledge, there are no incompatibilities between macOS, Microsoft Windows and Linux operating systems.

### Script 1: Reading files and building the dataset

The script `1RealData.R` reads multiple input datasets, including county-level AQI, economic indicators, population data and geographic characteristics, and merges them into a unified dataset suitable for our analysis. It handles missing values, standardizes identifiers, creates auxiliary variables and outputs two final CSV files: `data.csv` containing the target variable and the auxiliary variables, and `data_aux.csv` containing all the auxiliary information.

### Script 2: Calculating distances between counties

The script `2DistanceCounties.R` computes pairwise geographic distances between U.S. counties using the Haversine formula based on county centroids. It generates a normalized distance matrix, a spatial weight matrix with exponential decay and identifies the closest available counties for any missing subregion. The script also visualizes empirical semivariances over distance classes and maps county-level distance and weight information.  It outputs four final CSV files: `county_ALL_distance_matrix_haversine.csv`, `county_distance_matrix_haversine.csv`, `county_weights_matrix_haversine.csv` and \newline `missing_counties_dist.csv`.


### Script 3: Cross-validation for the spatio-temporal weights

The script `3STMQmodel_CV.R` performs leave-one-county-out cross-validation for STMQ models, combining spatial and temporal borrowing to predict median AQI values. It evaluates a grid of spatial decay parameters to select the optimal weighting. The auxiliary scripts for fitting, prediction and mean squared error (MSE) estimation (`QRLM.R` and `QRLMweights.R`) are used.


### Script 4: Fitting MQ, TMQ and STMQ models in the application to real data

The script `4STMQmodel_RealData.R` fits M-quantile (MQ), temporal MQ (TMQ) and spatio-temporal MQ (STMQ) models to county-level AQI data. It loads libraries and input files, standardizes covariates, and handles missing counties. The workflow fits MQ models, incorporates temporal and spatial borrowing for TMQ and STMQ models, and computes predictions for observed and out-of-sample counties and years. MSEs and bias-corrected predictions are calculated, and all results are output for further analysis and plotting. The auxiliary scripts for fitting, prediction and MSE estimation (`QRLM.R` and `QRLMweights.R`) are used.


### Script 5: Plotting and mapping

The script `5PlottingAndMapping.R` is used to visualize the results of the STMQ models. It produces plots comparing observed and predicted AQI values across counties and over time. It generates boxplots, scatterplots and time-series plots to assess model fit and trends. Spatial maps display predicted values and associated uncertainties, highlighting geographic patterns. The outputs are saved as high-quality PNGs for easy interpretation and reporting.
