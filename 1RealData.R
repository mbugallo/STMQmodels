
##################################################
###   Reading files and building the dataset   ###
##################################################

# INPUT DATASETS:
#  - Target_variables/annual_aqi_by_county_YYYY.csv : county-level AQI (target variable)
#  - Aux1_SAIPE/estYYall.csv                        : county-level poverty and income (SAIPE)
#  - Aux2_CAGDP/CAGDP1__ALL_AREAS_2001_2023.csv     : county-level GDP data
#  - Aux3_POP_State/ACSST1YYYY.S0101-Data.csv       : state-level population data
#  - Aux3_POP/ACSST1YYYY.S0101-Data.csv             : county-level population data
#  - Aux4_CoastlineCounties/coastline_counties.csv  : coastal county indicator
#  - Aux4_CoastlineCounties/Ruralurbancontinuumcodes2023.csv : metropolitan county indicator

# OUTPUT DATASETS:
#  - data.csv     : final analysis dataset including AQI and auxiliary variables
#  - data_aux.csv : auxiliary-only dataset (all covariates, no target variables)

# NOTE:
#  - Missing values in numerical variables are imputed using year-specific medians
#  - Some territories and non-continental states are excluded from the final datasets


# Load spatial and data manipulation libraries
suppressWarnings(
  suppressMessages({
    library(tigris)
    library(dplyr)
}))  

# Remove all objects from the current workspace to avoid conflicts
rm(list=ls())

# Create a lookup table linking state postal codes to full state names
postal_state  <- list(
  "AK" = "Alaska", 
  "AL" = "Alabama", 
  "AR" = "Arkansas", 
  "AZ" = "Arizona",
  "CA" = "California", 
  "CO" = "Colorado", 
  "CT" = "Connecticut",
  "DC" = "District of Columbia", 
  "DE" = "Delaware", 
  "FL" = "Florida",
  "GA" = "Georgia", 
  "HI" = "Hawaii", 
  "IA" = "Iowa", 
  "ID" = "Idaho",
  "IL" = "Illinois", 
  "IN" = "Indiana", 
  "KS" = "Kansas", 
  "KY" = "Kentucky",
  "LA" = "Louisiana", 
  "MA" = "Massachusetts", 
  "MD" = "Maryland",
  "ME" = "Maine", 
  "MI" = "Michigan", 
  "MN" = "Minnesota", 
  "MO" = "Missouri",
  "MS" = "Mississippi", 
  "MT" = "Montana", 
  "NC" = "North Carolina",
  "ND" = "North Dakota", 
  "NE" = "Nebraska", 
  "NH" = "New Hampshire",
  "NJ" = "New Jersey", 
  "NM" = "New Mexico", 
  "NV" = "Nevada",
  "NY" = "New York", 
  "OH" = "Ohio", 
  "OK" = "Oklahoma", 
  "OR" = "Oregon",
  "PA" = "Pennsylvania", 
  "PR" = "Puerto Rico", 
  "RI" = "Rhode Island",
  "SC" = "South Carolina", 
  "SD" = "South Dakota", 
  "TN" = "Tennessee",
  "TX" = "Texas", 
  "UT" = "Utah", 
  "VA" = "Virginia", 
  "VT" = "Vermont",
  "WA" = "Washington", 
  "WI" = "Wisconsin", 
  "WV" = "West Virginia",
  "WY" = "Wyoming")

postal_state <- data.frame(Postal.Code = names(postal_state), 
  State = unlist(postal_state), stringsAsFactors = FALSE)


# 1. TARGET VARIABLES

# Read all county-level AQI files and keep only the variables of interest
file_list <- list.files(path = "Target_variables", 
       pattern = "annual_aqi_by_county_.*\\.csv$",  full.names = TRUE)
data_target <- lapply(file_list, read.csv)
data_target <- lapply(data_target, function(df) {
  df[c('State', 'County', 'Year', 'Median.AQI')] }) 
  
# Attach state postal codes to the target data for consistent identifiers
data_target <- lapply(data_target, function(df) {
  merge(df, postal_state, by = "State", all.x = TRUE) })  


# 2. AUXILIARY VARIABLES I

# Read SAIPE auxiliary data (poverty and income) for all available years
file_list <- list.files(path = "Aux1_SAIPE", pattern = "est\\d{2}all\\.csv$", 
                        full.names = TRUE)
data_auxI <- lapply(file_list, read.csv, sep=';')

postal_state2  <- unique(data_auxI[[1]][, c('State.FIPS.Code', 'Postal.Code')])
names(postal_state2) <- c('State.Code', 'Postal.Code')
postal_state2$State.Code <- sprintf("%02d", postal_state2$State.Code)

# Standardize state and county FIPS codes and build a unique county identifier
# Remove state-level aggregates (code 000)
# Convert poverty percentages from character to numeric format
data_auxI  <- lapply(data_auxI, function(df) {	
	df[,1] <- sprintf("%02d", df[,1])
	df[,2] <- sprintf("%03d", df[,2])
	df <- df[df[, 2] != "000", ] 
	df$GeoFIPS <- paste0(df[,1], df[,2])
    df$County <- gsub(" County$", "", df$Name)  
    df <- df[c('GeoFIPS', 'County', 'Postal.Code', 
               'Poverty.Percent..All.Ages', 
               'Poverty.Percent..Age.0.17', 
               'Median.Household.Income')]	
    df$Poverty.Percent..All.Ages <- as.numeric(
      gsub(",", ".",  gsub("[^0-9,]", "", 
                           trimws(df$Poverty.Percent..All.Ages)))  )
    df$Poverty.Percent..Age.0.17 <- as.numeric(
      gsub(",", ".",  gsub("[^0-9,]", "", 
                           trimws(df$Poverty.Percent..Age.0.17)))  )
	return(df)	})
	
# Assign the corresponding year to each auxiliary dataset
years <- 2016:2023	
data_auxI <- lapply(1:length(data_auxI), function(i) {
  df <- data_auxI[[i]];  df$Year <- rep(years[i], dim(df)[1]);  return(df) 	})	


# 3. AUXILIARY VARIABLES II

data_auxII <- read.csv("Aux2_CAGDP/CAGDP1__ALL_AREAS_2001_2023.csv")
data_auxII$County <- sub(",.*", "", iconv(data_auxII$GeoName, 
                                          from = "latin1", to = "UTF-8"))
data_auxII2016 <- data_auxII[c('GeoFIPS', 'County', 'Description', 'LineCode', 'X2016')]
data_auxII2017 <- data_auxII[c('GeoFIPS', 'County', 'Description', 'LineCode', 'X2017')]
data_auxII2018 <- data_auxII[c('GeoFIPS', 'County', 'Description', 'LineCode', 'X2018')]
data_auxII2019 <- data_auxII[c('GeoFIPS', 'County', 'Description', 'LineCode', 'X2019')]
data_auxII2020 <- data_auxII[c('GeoFIPS', 'County', 'Description', 'LineCode', 'X2020')]
data_auxII2021 <- data_auxII[c('GeoFIPS', 'County', 'Description', 'LineCode', 'X2021')]
data_auxII2022 <- data_auxII[c('GeoFIPS', 'County', 'Description', 'LineCode', 'X2022')]
data_auxII2023 <- data_auxII[c('GeoFIPS', 'County', 'Description', 'LineCode', 'X2023')]

names(data_auxII2016) <- names(data_auxII2017) <- names(data_auxII2018) <- 
		names(data_auxII2019) <- names(data_auxII2020) <- names(data_auxII2021) <- 
		names(data_auxII2022) <- names(data_auxII2023) <- c('GeoFIPS', 'County', 
		'Description', 'LineCode', 'Variable')
data_auxII <- rbind(data_auxII2016, data_auxII2017, data_auxII2018, data_auxII2019, 
                    data_auxII2020, data_auxII2021, data_auxII2022, data_auxII2023)
data_auxII$Year <- c(rep('2016', dim(data_auxII2016)[1]), 
                     rep('2017', dim(data_auxII2017)[1]),
                     rep('2018', dim(data_auxII2018)[1]), 
                     rep('2019', dim(data_auxII2019)[1]),
                     rep('2020', dim(data_auxII2020)[1]), 
                     rep('2021', dim(data_auxII2021)[1]),
                     rep('2022', dim(data_auxII2022)[1]), 
                     rep('2023', dim(data_auxII2023)[1]))	 

# Reshape GDP data from wide to long format across years
data_auxII <- cbind(data_auxII[data_auxII$LineCode==1,], 
                    data_auxII$Variable[data_auxII$LineCode==2],
                    data_auxII$Variable[data_auxII$LineCode==3])[, c(1:2, 5:8)]
data_auxII$GeoFIPS <- trimws(data_auxII$GeoFIPS)
names(data_auxII) <- c('GeoFIPS', 'County', 'Real.GDP', 'Year', 
                       'Quantity.Index', 'Current.Dollar.GDP')
 
data_auxII <- data_auxII[substr(as.character(data_auxII$GeoFIPS), 3, 5) != "000", ]
data_auxII$State.Code <- substr(as.character(data_auxII$GeoFIPS), 1, 2) 

data_auxII <- merge(data_auxII, postal_state2, by = "State.Code", all.x = TRUE)

common_GeoFIPS <- unique(Reduce(union, c(list(data_auxII$GeoFIPS), 
                         lapply(data_auxI, function(df) df$GeoFIPS))))
common_GeoFIPS <- common_GeoFIPS[!is.na(common_GeoFIPS)]
common_counties <- unique(Reduce(union, c(list(data_auxII$County), 
                          lapply(data_auxI, function(df) df$County))))
common_counties <- common_counties[!is.na(common_counties)]

# Convert GDP variables to numeric and scale units
data_auxII$Real.GDP <- as.numeric(gsub("[^0-9.]", "", 
                                 trimws(data_auxII$Real.GDP)))/1000
data_auxII$Current.Dollar.GDP <- as.numeric( gsub("[^0-9.]", "", 
                                 trimws(data_auxII$Current.Dollar.GDP)))/1000

# 4. AUXILIARY VARIABLES III - State level variables

# Read state-level population data and compute population densities
# 2020 is missing. Create Aux3_POP_State/ACSST1Y2020.S0101-Data.csv equal to that of 2019
file_list <- list.files(path = "Aux3_POP_State", 
                        pattern = "ACSST1Y20[0-9]{2}\\.S0101-Data", 
                        full.names = TRUE)
data_auxIII <- lapply(file_list, read.csv, sep=';')
data_auxIII <- lapply(data_auxIII, function(df) {
  df <- df[, 1:3, drop = FALSE] 
  colnames(df) <- c("State.Code", "State", "Total.Pop.State")  
  df$State.Code <- sub("0400000US", "", df$State.Code); df})

data_auxIII <- lapply(1:length(data_auxIII), function(i) {
  df <- data_auxIII[[i]];  df$Year <- rep(years[i], dim(df)[1]);  return(df) 	})	

data_auxIII <-  do.call(rbind, data_auxIII)

states <- data.frame(states(cb = TRUE))[, c("GEOID", "ALAND")]
states$ALAND <- states$ALAND/1e6 # resultados en m2
names(states) <- c("State.Code", 'Area.Land.State')
  
# Merge land area information to compute state population densities
data_auxIII <- merge(data_auxIII, states, by="State.Code", all.x=T)
data_auxIII$Dens.Pop.State <- 
  data_auxIII$Total.Pop.State/data_auxIII$Area.Land.State

# 4. AUXILIARY VARIABLES III - County level variables

# Read county-level population data and align it with the common county set
file_list <- list.files(path = "Aux3_POP", 
                        pattern = "ACSST1Y20[0-9]{2}\\.S0101-Data", 
                        full.names = TRUE)
data_auxIII.C <- lapply(file_list, read.csv, sep=';')
data_auxIII.C <- lapply(data_auxIII.C, function(df) {
  df <- df[, c(1,3), drop = FALSE] 
  colnames(df) <- c("County.Code", "Total.Pop.County")  
  df$County.Code <- sub("0500000US", "", df$County.Code); df})

data_auxIII.C <- lapply(data_auxIII.C, function(df) 
  merge(data.frame('County.Code' = common_GeoFIPS), 
        df, by='County.Code', all.x=T))

data_auxIII.C <- lapply(1:length(data_auxIII.C), function(i) {
  df <- data_auxIII.C[[i]] 
  df$Year <- rep(years[i], dim(df)[1]);  return(df) 	})

data_auxIII.C <-  do.call(rbind, data_auxIII.C)

data_auxIII.C$State.Code <- substr(data_auxIII.C$County.Code, 1, 2)
data_auxIII.C <- merge(data_auxIII.C, data_auxIII, 
                       by=c('State.Code', 'Year'), all.x=T)

# Impute missing county populations by distributing the state residual equally
sum.NAN.county <- data.frame(data_auxIII.C %>% group_by(State.Code, Year) %>%
         summarise(NA.Count = sum(is.na(Total.Pop.County)), .groups = "drop"))

sum.pop.county <- aggregate(Total.Pop.County ~ State.Code + Year, 
                            data = data_auxIII.C, FUN = sum, na.rm=T)
names(sum.pop.county) <- c('State.Code', 'Year', 'Total.Pop.County.Sample')
sum.pop.county <- merge(sum.pop.county, sum.NAN.county, by=c('State.Code', 'Year'))
sum.pop.county <- merge(data_auxIII, sum.pop.county,
                        by=c('State.Code', 'Year'), all.x=T, all.y=T)
sum.pop.county$sharing <- round((sum.pop.county$Total.Pop.State-
                                   sum.pop.county$Total.Pop.County)/sum.pop.county$NA.Count,0)
sum.pop.county$sharing[is.na(sum.pop.county$sharing )] <- 0

data_auxIII.C <- merge(data_auxIII.C, sum.pop.county[, c('State.Code', 'Year', 'sharing')],
                       by=c('State.Code', 'Year'), all.x=T)
data_auxIII.C$Total.Pop.County[is.na(data_auxIII.C$Total.Pop.County)] <- 
  data_auxIII.C$sharing[is.na(data_auxIII.C$Total.Pop.County)]

counties <- data.frame(counties(cb = TRUE))[, c("GEOID", "ALAND")]
counties$ALAND <- counties$ALAND/1e6 # resultados en m2
names(counties) <- c("County.Code", 'Area.Land.County')

data_auxIII.C <- merge(data_auxIII.C, counties, by="County.Code", all.x=T)
data_auxIII.C$Dens.Pop.County <- data_auxIII.C$Total.Pop.County/data_auxIII.C$Area.Land.County

data_auxIII.C <- data_auxIII.C[, c('County.Code', 'Year', 
                                   'Total.Pop.County', 'Dens.Pop.County')]
names(data_auxIII.C)[1] <- 'GeoFIPS'

# 5. AUXILIARY VARIABLES IV

data_auxIV <- read.csv('Aux4_CoastlineCounties/coastline_counties.csv', sep=';')[,1]
# Create a binary indicator for coastal counties
data_auxIV <- data.frame('GeoFIPS'= sprintf("%05d", data_auxIV), 'ID.Coast' = 1)


# 6. AUXILIARY VARIABLES V

data_auxV <- read.csv('Aux4_CoastlineCounties/Ruralurbancontinuumcodes2023.csv', sep=',')
data_auxV <- data_auxV[data_auxV$Attribute=='RUCC_2023', ]
# Create a metropolitan indicator based on the RUCC classification
data_auxV$ID.Metrop <- (data_auxV$Value=='1'|data_auxV$Value=='2'|data_auxV$Value=='3')*1
 
data_auxV <- data_auxV[, c('FIPS', 'ID.Metrop')]
data_auxV$FIPS <- sprintf("%05d", data_auxV$FIPS)
names(data_auxV) <- c('GeoFIPS', 'ID.Metrop')


# 7. COMMON COUNTIES
# Keep only counties common to target and auxiliary datasets
data_target <- Reduce('rbind', lapply(data_target, 
                      function(df) df[df$County %in% common_counties, ]))
# Build a unique identifier combining state, county and year
data_target$Identifier <- paste(data_target$Postal.Code, 
                                data_target$County, data_target$Year, sep = "_")

data_auxI <- Reduce('rbind', lapply(data_auxI, 
                             function(df) df[df$County %in% common_counties, ]))
data_auxI$Identifier <- paste(data_auxI$Postal.Code, 
                              data_auxI$County, data_auxI$Year, sep = "_")

data_auxII <- data_auxII[data_auxII$County %in% common_counties, ]
data_auxII$Identifier <- paste(data_auxII$Postal.Code, 
                               data_auxII$County, data_auxII$Year, sep = "_")

County.id  <- unique(data_auxII[, c('GeoFIPS', 'County')])
GeoFIPS.id <- County.id[County.id$County %in% common_counties, 1]


# 8. COMMON TARGET AND AUXILIARY VARIABLES: output data.csv 
# Merge target and auxiliary variables into a single final dataset
data_auxI_II <- merge(data_auxI[, c('GeoFIPS', 'Poverty.Percent..All.Ages', 'Year',
                'Poverty.Percent..Age.0.17', 'Median.Household.Income')],
                data_auxII[, c('GeoFIPS', 'State.Code', 'Real.GDP', 
                'Current.Dollar.GDP', 'Identifier', 'Year')],
                by=c('GeoFIPS', 'Year'), all.x=TRUE)

final.set <- merge(merge(merge(merge(data_target, data_auxI_II[, c('GeoFIPS', 
            'Poverty.Percent..All.Ages', 'Poverty.Percent..Age.0.17', 
            'Median.Household.Income', 'Identifier', 'State.Code', 'Real.GDP', 
            'Current.Dollar.GDP')], by='Identifier', all.x=T), 
            data_auxIII[, c('State.Code', 'Total.Pop.State', 
                            'Dens.Pop.State', 'Year')], 
            by=c('State.Code', 'Year'), all.x=T), 
            data_auxIII.C, by=c('GeoFIPS', 'Year')),
            data_auxV, by='GeoFIPS', all.x=T)

final.set <- final.set[, c('GeoFIPS', 'County', 'State.Code', 'Postal.Code', 
        'State', 'Year', 'Median.AQI', 'Median.Household.Income',
        'Poverty.Percent..All.Ages', 'Poverty.Percent..Age.0.17', 
        'Real.GDP', 'Current.Dollar.GDP', 'Total.Pop.County', 'Dens.Pop.County',
        'Total.Pop.State', 'Dens.Pop.State', 'ID.Metrop')]
				
final.set[, c('Median.AQI', 'Median.Household.Income', 'Poverty.Percent..All.Ages',
		'Poverty.Percent..Age.0.17', 'Real.GDP', 'Current.Dollar.GDP',
		'Total.Pop.County', 'Dens.Pop.County', 'Total.Pop.State', 'Dens.Pop.State')] <-
		apply(final.set[, c('Median.AQI', 'Median.Household.Income', 
		                 'Poverty.Percent..All.Ages', 'Poverty.Percent..Age.0.17', 
		                 'Real.GDP', 'Current.Dollar.GDP',	'Total.Pop.County',
		                 'Dens.Pop.County', 'Total.Pop.State', 'Dens.Pop.State')],
		      2, as.numeric)

# Impute remaining missing values using year-specific medians
final.set <- final.set %>% group_by(Year) %>% mutate(across(7:16, ~ {
  med <- median(.x, na.rm = TRUE); .x[is.na(.x)] <- med; return(.x) }))

final.set <- merge(final.set, data_auxIV, by = 'GeoFIPS', all.x = TRUE)
final.set$ID.Coast[is.na(final.set$ID.Coast)] <- 0

# Exclude non-continental states and territories
final.set <- final.set[apply(final.set, 1, function(x) all(!is.na(x))), ]
final.set <- final.set %>% arrange(Year, GeoFIPS)

final.set <- final.set[!final.set$Postal.Code%in%c('AK', 'AS', 'MP', 
                                  'DC', 'PR', 'RI', 'HI', 'VI', 'GU'),]
final.set <- final.set %>% arrange(Year, GeoFIPS)

# AK- Alaska (state)
# AS - American Samoa (land), 
# MP - Northern Mariana Islands (land)
# PR - Puerto Rico (land)
# RI - Rhode Island (state)
# HI - Hawaii (state)
# VI - U.S. Virgin Islands (land) 
# GU - Guam (land)
# DC - District of Columbia (land)

# Convert poverty rates from percentages to proportions
final.set[, c('Poverty.Percent..All.Ages', 'Poverty.Percent..Age.0.17')] <-
  final.set[, c('Poverty.Percent..All.Ages', 'Poverty.Percent..Age.0.17')]/100
names(final.set)[names(final.set) %in% c('Poverty.Percent..All.Ages', 
             'Poverty.Percent..Age.0.17')] <- 
  c('Poverty.Proportion.All.Ages', 'Poverty.Proportion.Age.0.17')

# Save the target data
write.csv(final.set, "data.csv", row.names = FALSE)


# 9. COMMON AUXILIARY VARIABLES: output data_aux.csv

final.set.aux <- merge(merge(merge(merge(data_auxI[, c('GeoFIPS', 'County', 'Postal.Code', 
                 'Year', 'Poverty.Percent..All.Ages', 'Poverty.Percent..Age.0.17', 
                 'Median.Household.Income', 'Identifier' )],  
                 data_auxII[, c('GeoFIPS', 'Year','State.Code', 'Real.GDP',
                                'Current.Dollar.GDP', 'Identifier')], 
                 by=c('GeoFIPS', 'Year'), all.x=T), data_auxIII[, c('State', 
                 'State.Code', 'Total.Pop.State', 'Dens.Pop.State', 'Year')], 
                 by=c('State.Code', 'Year'), all.x=T), 
                 data_auxIII.C, by=c('GeoFIPS', 'Year')),
                 data_auxV, by='GeoFIPS', all.x=T)

final.set.aux <- final.set.aux[, c('GeoFIPS', 'County', 'State.Code', 'Postal.Code', 
                'State', 'Year', 'Median.Household.Income', 'Poverty.Percent..All.Ages', 
                'Poverty.Percent..Age.0.17', 'Real.GDP', 'Current.Dollar.GDP', 
                'Total.Pop.County', 'Dens.Pop.County', 'Total.Pop.State', 
                'Dens.Pop.State', 'ID.Metrop')]

clean_numeric <- function(x) { as.numeric(gsub("[^0-9.-]", "", trimws(x))) }

final.set.aux[, c('Median.Household.Income', 'Poverty.Percent..All.Ages', 
                  'Poverty.Percent..Age.0.17', 'Real.GDP', 'Current.Dollar.GDP', 
                  'Total.Pop.County', 'Dens.Pop.County', 'Total.Pop.State', 'Dens.Pop.State')] <-
  apply(final.set.aux[, c('Median.Household.Income', 'Poverty.Percent..All.Ages', 
                  'Poverty.Percent..Age.0.17', 'Real.GDP',  'Current.Dollar.GDP', 
                  'Total.Pop.County', 'Dens.Pop.County', 'Total.Pop.State', 
                  'Dens.Pop.State')], 2, clean_numeric)

final.set.aux <- final.set.aux %>% group_by(Year) %>% mutate(across(7:15, ~ {
  med <- median(.x, na.rm = TRUE); .x[is.na(.x)] <- med; return(.x) }))

final.set.aux <- merge(final.set.aux, data_auxIV, by = 'GeoFIPS', all.x = TRUE)
final.set.aux$ID.Coast[is.na(final.set.aux$ID.Coast)] <- 0

final.set.aux <- final.set.aux[apply(final.set.aux, 1, function(x) all(!is.na(x))), ]
final.set.aux <- final.set.aux %>% arrange(Year, GeoFIPS)

final.set.aux <- final.set.aux[!final.set.aux$Postal.Code%in%c('AK', 'AS', 'MP', 
                                          'DC', 'PR', 'RI', 'HI', 'VI', 'GU'),]
final.set.aux <- final.set.aux %>% arrange(Year, GeoFIPS)

final.set.aux[, c('Poverty.Percent..All.Ages', 'Poverty.Percent..Age.0.17')] <-
  final.set.aux[, c('Poverty.Percent..All.Ages', 'Poverty.Percent..Age.0.17')]/100
names(final.set.aux)[names(final.set.aux) %in% c('Poverty.Percent..All.Ages', 
            'Poverty.Percent..Age.0.17')] <-  c('Poverty.Proportion.All.Ages', 
                                                'Poverty.Proportion.Age.0.17')

# Save the auxiliary information
write.csv(final.set.aux, "data_aux.csv", row.names = FALSE)

