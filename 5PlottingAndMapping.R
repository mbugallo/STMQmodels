
#####################################################
###   Plotting results and mapping predictions    ###
#####################################################

# Remove all objects from the current workspace to avoid conflicts
rm(list=ls())

# Load the workspace
load('STMQpred.Rdata')

# Load required libraries
suppressWarnings(
  suppressMessages({
    attach(data)
    library(dplyr)
    library(reshape2)
    library(ggplot2)
    library(sf)
    library(tigris)
    library(ggspatial)
    library(viridis)
    library(tidyr)
}))

options(tigris_use_cache = TRUE)

color.legend <- c("0-5" = "#08306b", "5-20" = "#2879a3", "20-40" = "#43a2ca", 
                  "40-60" = "#7bccc4", "60-80" = "#c2e699", "80-95" = "#fecc5c", 
                  "95-100" = "#e31a1c", "NA" = "grey")

probs <- c(0, 0.05, 0.2, 0.4, 0.6, 0.8, 0.95, 1)
labels <- c("0-5", "5-20", "20-40", "40-60",  "60-80", "80-95", "95-100")


##########################################
###    Spaghetti plot with ggplot2     ###
##########################################

dt_year_quantiles <- pred.rt.id %>% group_by(Year) %>%
  summarise(Q05 = quantile(pred.rt, 0.05, na.rm = TRUE),
            Q25 = quantile(pred.rt, 0.25, na.rm = TRUE),
            Q50 = quantile(pred.rt, 0.50, na.rm = TRUE),
            Q75 = quantile(pred.rt, 0.75, na.rm = TRUE),
            Q95 = quantile(pred.rt, 0.95, na.rm = TRUE))

dt_state <- pred.rt.id %>% group_by(State.Code, Year) %>%
  summarise(pred_med = median(pred.rt, na.rm = TRUE), .groups = "drop")

ggplot(dt_state, aes(Year, pred_med, group = State.Code)) +
  geom_line(alpha = 0.85, colour = "grey50") +
  stat_summary(aes(group = 1), fun = median, geom = "line",
    linewidth = 1.3, col='red') +
  labs(x = "Year", y = "Median predicted AQI") + theme_minimal(base_size = 18) +
  theme(axis.title.x = element_text(size = 18, margin = margin(t = 18)),
    axis.title.y = element_text(size = 18, margin = margin(r = 18)),
    axis.text   = element_text(size = 18),
    panel.grid.minor = element_blank()
  )

# ggsave('AQIspag.jpg', quality = 80, width = 6, height = 6, units = 'in', dpi = 150)

##################################################
###       Plotting MQ - Ranks (STATES)       #####
###         OUT-OF-SAMPLE COUNTIES           #####
##################################################

# Load U.S. states
states <- states(cb = TRUE)
year.plot <- '2016'

# Merge the shapefile with the GeoFIPS data
pred.color.r <- dt_state[dt_state$Year==year.plot, ]

states <- states %>% left_join(pred.color.r, by = c("GEOID" = "State.Code"))

states_sf <- st_as_sf(states)
states_sf$pred_category <- cut(states_sf$pred_med, breaks = quantile(states_sf$pred_med, 
                               probs = probs,  na.rm = TRUE), labels = labels, 
                                 include.lowest = TRUE)

quartz(width = 9, height = 6) 
ggplot(data = states_sf) +
  geom_sf(aes(fill = pred_category), color = NA, size = 0) +
  scale_fill_manual(name = "Percentile", values = color.legend,
                    na.value = "grey",  na.translate = FALSE) +
  theme_minimal() +
  labs(  title = "", fill = "Percentile") +
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

# ggsave('AQIstaterdt.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)

################################
###    Plotting Counties     ###
################################

# Load U.S. counties
counties <- counties(cb = TRUE) 

# Merge the shapefile with the GeoFIPS data
GeoFIPS.id <- data.frame('GEOID'=sprintf("%05d", GeoFIPS), 
                         'color'=rep('red', length(GeoFIPS)))
counties <- counties %>% left_join(GeoFIPS.id, by = c("GEOID" = "GEOID"))

# Convert to sf object if needed
counties_sf <- st_as_sf(counties)

# Plot the map with colors based on the GeoFIPS data
quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = factor(color)), color = "white", size = 0) + 
  theme_minimal() +
  labs(title = " ") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

# ggsave('Counties.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)


#######################################
###       Plotting MQ - Counties    ###
#######################################

counties <- counties(cb = TRUE) 
year.plot <- '2023'

# Merge the shapefile with the GeoFIPS data
pred.color.s <- data.frame('GEOID'=sprintf("%05d", 
                            pred.st.id$GeoFIPS)[pred.st.id$Year==year.plot], 
                           'pred.st'=pred.st.id$pred.st[pred.st.id$Year==year.plot])


counties <- counties %>% left_join(pred.color.s, by = c("GEOID" = "GEOID"))
counties_sf <- st_as_sf(counties)

# Plot the map with colors based on the GeoFIPS data
quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred.st), color = "white", size = 0.1) + 
  scale_fill_viridis_c(option = "plasma", name = "AQI Scale") + 
  theme_minimal() +
  labs(title = " ", fill = "Variable") +
  theme(legend.position = c(0.9, 0.25)) +
  xlim(-125, -65) + ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)


#####################################
###      Plotting MQ - Ranks    #####
###    OUT-OF-SAMPLE COUNTIES   #####
#####################################

counties <- counties(cb = TRUE)  
year.plot <- '2023'

# Merge the shapefile with the GeoFIPS data
pred.color.r <- data.frame('GEOID'=sprintf("%05d", 
                            pred.rt.id$GeoFIPS)[pred.rt.id$Year==year.plot], 
                           'pred.rt'=pred.rt.id$pred.rt[pred.rt.id$Year==year.plot])

counties <- counties %>% left_join(pred.color.r, by = c("GEOID" = "GEOID"))
counties_sf <- st_as_sf(counties)

counties_sf$pred_category <- cut(counties_sf$pred.rt, breaks = quantile(counties_sf$pred.rt, 
                                 probs = probs,  na.rm = TRUE), labels = labels, 
                                 include.lowest = TRUE)

quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred_category), color = NA, size = 0) +
  scale_fill_manual(name = "Percentile", values = color.legend, na.value = "grey" ) +
  theme_minimal() +
  labs(  title = "", fill = "Percentile") +
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

# ggsave('AQIrdt.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)


#####################################
###       Plotting MQ - Ranks   #####
###         SAMPLE COUNTIES     #####
#####################################

counties_sf.drop <- st_drop_geometry(counties_sf)
counties_sf <- st_as_sf(merge(merge(counties_sf.drop, pred.color.s,  by='GEOID', 
                              sort=FALSE), data.frame(counties_sf[, c('GEOID', 'geometry')]), 
                              by='GEOID', sort=FALSE, all.y=TRUE))

quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred_category), color = NA, size = 0) +
  scale_fill_manual( name = "Percentile", values = color.legend, na.value = "grey") +
  theme_minimal() +
  labs(title = " ",
       subtitle = " ",
       fill = "Percentile") +    
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

# ggsave('AQIsdt.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)


################################################
###       Plotting TEMPORAL MQ - Counties    ###
################################################

counties <- counties(cb = TRUE) 
year.plot <- '2023'

# Merge the shapefile with the GeoFIPS data
pred.color.s <- data.frame('GEOID'=sprintf("%05d", 
                                           pred.TMQ.st.id$GeoFIPS)[pred.TMQ.st.id$Year==year.plot], 
                           'pred.st'=pred.TMQ.st.id$pred.st[pred.TMQ.st.id$Year==year.plot])


counties <- counties %>% left_join(pred.color.s, by = c("GEOID" = "GEOID"))
counties_sf <- st_as_sf(counties)

# Plot the map with colors based on the GeoFIPS data
quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred.st), color = "white", size = 0.1) + 
  scale_fill_viridis_c(option = "plasma", name = "AQI Scale") + 
  theme_minimal() +
  labs(title = " ", fill = "Variable") +
  theme(legend.position = c(0.9, 0.25)) +
  xlim(-125, -65) + ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)


##############################################
###      Plotting TEMPORAL MQ - Ranks    #####
###         OUT-OF-SAMPLE COUNTIES       #####
##############################################

counties <- counties(cb = TRUE)  
year.plot <- '2023'

# Merge the shapefile with the GeoFIPS data
pred.color.r <- data.frame('GEOID'=sprintf("%05d", 
                            pred.TMQ.rt.id$GeoFIPS)[pred.TMQ.rt.id$Year==year.plot], 
                           'pred.rt'=pred.TMQ.rt.id$pred.rt[pred.TMQ.rt.id$Year==year.plot])

counties <- counties %>% left_join(pred.color.r, by = c("GEOID" = "GEOID"))
counties_sf <- st_as_sf(counties)

counties_sf$pred_category <- cut(counties_sf$pred.rt, 
                                 breaks = quantile(counties_sf$pred.rt, 
                                 probs = probs,  na.rm = TRUE), labels = labels, 
                                 include.lowest = TRUE)

quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred_category), color = NA, size = 0) +
  scale_fill_manual(name = "Percentile", values = color.legend, na.value = "grey" ) +
  theme_minimal() +
  labs(  title = "", fill = "Percentile") +
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

# ggsave('AQIrdtTMQ.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)


##############################################
###       Plotting TEMPORAL MQ - Ranks   #####
###           SAMPLE COUNTIES            #####
##############################################

counties_sf.drop <- st_drop_geometry(counties_sf)
counties_sf <- st_as_sf(merge(merge(counties_sf.drop, pred.color.s,  by='GEOID', 
                                    sort=FALSE), data.frame(counties_sf[, c('GEOID', 'geometry')]), 
                              by='GEOID', sort=FALSE, all.y=TRUE))

quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred_category), color = NA, size = 0) +
  scale_fill_manual( name = "Percentile", values = color.legend, na.value = "grey") +
  theme_minimal() +
  labs(title = " ",
       subtitle = " ",
       fill = "Percentile") +    
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

# ggsave('AQIsdtTMQ.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)


###############################################
###    Boxplots BETA's Temporal GWR MQ      ###
###############################################

data_long <- as.data.frame(beta.SAE.STMQ) %>%
  pivot_longer(cols = 3:(p+2), names_to = "Variable", values_to = "Value") %>%
  mutate(
    Variable = factor(Variable, levels = names(as.data.frame(beta.SAE.STMQ)[3:(p+2)])),
    Year = factor(Year)
  )

quartz(width = 10, height = 6)
ggplot(data_long, aes(x = Year, y = Value, fill = Variable)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  
  # geom_hline(yintercept = 0, color = "black", linetype = "solid") + 
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c(2016, '', 2018, '', 2020, '', 2022, '')) + 
  labs(x = "Year", y = "Estimation") + 
  theme_minimal() +  
  theme(
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(size = 14),
    panel.spacing = unit(0.5, "lines")
  ) + 
  facet_wrap(~Variable, ncol = 4, scales = "free_y") + 
  guides(fill = guide_legend(nrow = 1))


#################################
###     Plotting weights      ###
#################################

df <- data.frame(
  Value = weights.df.counties[[154]] * ww[[1]],  
  Index = seq_along(weights.df.counties[[154]] * ww[[1]])
)

# Crear el gráfico
quartz(width = 6, height = 6)
ggplot(df, aes(x = Index, y = Value)) + 
  geom_point(alpha = 0.3, color = "black", size=2) + 
  labs(x = "Year", y = expression(w[dtj * "," * gil])) + 
  theme_minimal() + 
  theme(
    axis.title.x = element_text(size = 25),  
    axis.title.y = element_text(size = 25),  
    axis.text.x = element_text(angle = 45, hjust = 1, size = 25),
    axis.text.y = element_text(size = 25),
    legend.position = "none",
    strip.text = element_text(size = 25),
    panel.spacing = unit(0.5, "lines")
  ) +  ylim(0, 1) + 
  scale_x_continuous(breaks = cumsum(as.numeric(table(Year))) - as.numeric(table(Year))[1]/2, 
                     labels = years) +  
  geom_vline(xintercept = cumsum(as.numeric(table(Year))) - as.numeric(table(Year))[1], 
             linetype = "dotted", color = "gray90") +  
  geom_hline(yintercept = pretty(weights.df.counties[[154]] * ww[[1]]), 
             linetype = "dotted", color = "gray90")



###########################################
###     Validation Temporal GWR MQ      ###
###########################################

quartz(width = 7, height = 6)
ggplot() +
  geom_boxplot(aes(x = pred.GWR$Postal.Code, y = res.st.GWR), 
               fill = "gray70", color = "darkblue") +
  geom_hline(yintercept = 3, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -3, color = "red", linetype = "dashed") +
  labs(title = " ", 
       x = "State", y = "Standardized residuals") +
  scale_x_discrete(breaks = function(x) quantile(1:m)) +  #
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18, hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16))

quartz(width = 7, height = 6)
ggplot() +
  geom_boxplot(aes(x = as.character(pred.GWR$Year), y = res.st.GWR), 
               fill = "gray70", color = "darkblue") +
  geom_hline(yintercept = 3, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -3, color = "red", linetype = "dashed") +
  labs(title = " ", 
       x = "Year", y = "Standardized residuals") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18, hjust = 1), 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16))

###############################################
###     Comparing MQ model predictions      ###
###############################################

pred.MQ.TMQ.STMQ <- merge(merge(pred.st.id[, c('Year', 'GEOID', 'id', 
                          'Postal.Code', 'Median.AQI', 'pred.st')], 
                          pred.TMQ.st.id[,c('id', 'pred.st')], by='id'), 
                          pred.GWR[,c('id', 'pred.st')], by='id')

names(pred.MQ.TMQ.STMQ) <- c('id', 'Year', 'GEOID', 'Postal.Code',
                             'AQI', 'MQ', 'TWMQ', 'STMQ')

dflong <- melt(pred.MQ.TMQ.STMQ, id.vars = "Year", 
               measure.vars = c('AQI', "MQ", "TWMQ", "STMQ"))

quartz(width = 11.5, height = 5)
ggplot(dflong, aes(x = factor(Year), y = value, fill = variable)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2, width = 0.6) +
  scale_fill_manual(values = c("#FFD700", "#1b9e77", "#d95f02", "#7570b3")) +  
  labs(title = " ",
       x = "Year",
       y = "Predictions",
       fill = "") +
  ylim(0,90)+
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 18), 
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 18, hjust = 1), 
    axis.title.y = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)
  )


###############################################
###     Comparing STMQ model predictions    ###
###############################################


pred.STMQs <- merge(data[, c('Year', 'GEOID', 'Median.AQI')], 
                    MSE[, c('Year', 'GEOID', 'Postal.Code', 
                            'predN.mse', 'predBC.mse')],
                    by=c('Year', 'GEOID'))
names(pred.STMQs) <- c('Year', 'GEOID', 'AQI', 'Postal.Code', 'STMQ', 'BSTMQ')

dflong <- melt(pred.STMQs, id.vars = "Year", 
               measure.vars = c('AQI', "STMQ", "BSTMQ"))

quartz(width = 11.5, height = 5)
ggplot(dflong, aes(x = factor(Year), y = value, fill = variable)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2, width = 0.6) +
  scale_fill_manual(values = c("#FFD700", "#7570b3", "lightcoral")) +  
  labs(title = " ",
       x = "Year",
       y = "Predictions",
       fill = "") +
  ylim(0,90)+
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 18), 
    legend.text = element_text(size = 18),
    axis.text.x = element_text(size = 18, hjust = 1), 
    axis.title.y = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)
  )

################################################
###    Plotting Temporal GWR MQ - Counties   ###
################################################

counties <- counties(cb = TRUE) 
year.plot <- '2023'

# Merge the shapefile with the GeoFIPS data
pred.color.s.GWR <- data.frame('GEOID'=sprintf("%05d", pred.GWR$GeoFIPS)[pred.GWR$Year==year.plot], 
                               'pred.st'=pred.GWR$pred.st[pred.GWR$Year==year.plot])


counties <- counties %>% left_join(pred.color.s.GWR, by = c("GEOID" = "GEOID"))
counties_sf <- st_as_sf(counties)

# Plot the map with colors based on the GeoFIPS data
quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred.st), color = "white", size = 0.1) + 
  scale_fill_viridis_c(option = "plasma", name = "AQI Scale") + theme_minimal() +
  labs(title = " ", fill = "Variable") +
  theme(legend.position = c(0.9, 0.25)) +
  xlim(-125, -65) + ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)


##############################################
###   Plotting Temporal GWR MQ - Ranks   #####
###       OUT-OF-SAMPLE COUNTIES         #####
##############################################

counties <- counties(cb = TRUE)  
year.plot <- '2023'

# Merge the shapefile with the GeoFIPS data
pred.color.r.GWR <- data.frame('GEOID'=sprintf("%05d", predr.GWR$GeoFIPS)[predr.GWR$Year==year.plot], 
                               'pred.rt'=predr.GWR$pred.rt[predr.GWR$Year==year.plot])


counties <- counties %>% left_join(pred.color.r.GWR, by = c("GEOID" = "GEOID"))
counties_sf <- st_as_sf(counties)

counties_sf$pred_category <- cut(counties_sf$pred.rt,
                                 breaks = quantile(counties_sf$pred.rt,  probs = probs, 
                                 na.rm = TRUE), labels = labels, include.lowest = TRUE)

quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred_category), color = NA, size = 0) +
  scale_fill_manual(name = "Percentile", values = color.legend, na.value = "grey") +
  theme_minimal() +
  labs(  title = "", fill = "Percentile") +
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +   ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

# ggsave('AQIrdtSTMQ.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)


##############################################
###   Plotting Temporal GWR MQ - Ranks   #####
###           SAMPLE COUNTIES            #####
##############################################

counties_sf.drop <- st_drop_geometry(counties_sf)
counties_sf <- st_as_sf(merge(merge(counties_sf.drop, pred.color.s, by='GEOID', 
                           sort=FALSE), data.frame(counties_sf[, c('GEOID', 'geometry')]), 
                           by='GEOID', sort=FALSE, all.y=TRUE))

quartz(width = 9, height = 6) 
ggplot(data = counties_sf) +
  geom_sf(aes(fill = pred_category), color = NA, size = 0) +
  scale_fill_manual( name = "Percentile", values = color.legend, na.value = "grey") +
  theme_minimal() +
  labs(title = " ",
       subtitle = " ",
       fill = "Percentile") +    
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 16),            
        axis.title.y = element_text(size = 16),            
        axis.text.x = element_text(size = 14),             
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),           
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, hjust = 0.5)) +    
  xlim(-125, -65) +    
  ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

# ggsave('AQIsdtSTMQ.jpg', quality = 80, width = 9, height = 6, units = 'in', dpi = 150)


########################################
###   Some results for the paper   #####
########################################

mqo.GEOID.STMQ <- data.frame(mqo.GEOID.STMQ, 
                  'State.Code'= substr(mqo.GEOID.STMQ$GEOID, 1, 2))

mqo.GEOID.STMQ <- merge(merge(mqo.GEOID.STMQ, unique(data[, c('Postal.Code', 
                        'State.Code')]), by='State.Code'), 
                        aggregate(as.numeric(qmat[,1]),by=list(qmat[,2]),mean),
                        by.x = 'Postal.Code', by.y = 'Group.1')

summary(100*(mqo.GEOID.STMQ[,4]- mqo.GEOID.STMQ[,5])/mqo.GEOID.STMQ[,5])


#############################################
###    MEAN SQUARED ERROR ESTIMATION    #####
#############################################

# FIRST PLOT

quartz(width = 7, height = 6) 
ggplot(MSE, aes(x = as.character(Year), y = 100*biasBC/predBC.mse)) +
  geom_boxplot(outlier.shape = 16, fill = "gray70", color = "darkblue") +
  labs(x = "Year", y = "RBIAS (%)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 18, hjust = 1),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))

# SECOND PLOT

quartz(width = 7, height = 6) 
ggplot(MSE, aes(x = as.character(Year), y = 100*rrmse3BC)) +
  geom_boxplot(outlier.shape = 16, fill = "gray70", color = "darkblue") +
  labs(x = "Year", y = "RRMSE (%)", title = "") +
  theme_minimal() + ylim(0, 45) +
  theme(axis.text.x = element_text(size = 18, hjust = 1),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))



################################################
###     Mapping RRMSE - Temporal GWR MQ      ###
################################################

counties <- counties(cb = TRUE) 
year.plot <- '2023'

# Merge the shapefile with the GeoFIPS data
pred.color.s.GWR <- data.frame('GEOID'=MSE$GEOID, 
                    'RRMSE'=100*MSE$rrmse1BC)[MSE$Year==year.plot, ] 
                    #& MSE$id.in.s==1


counties <- counties %>% left_join(pred.color.s.GWR, by = c("GEOID" = "GEOID"))
counties_sf <- st_as_sf(counties)

counties_sf$RRMSE <- cut(counties_sf$RRMSE,  
                         breaks = c(0, 5, 10, 20, 30, 40, Inf),  
                         labels = c("0-5", "5-10", "10-20", "20-30", "30-40", "40+"),  
                         include.lowest = TRUE)  

quartz(width = 9, height = 6)  
ggplot(data = counties_sf) +  
  geom_sf(aes(fill = RRMSE), color = "white", size = 0.1) +  
  scale_fill_manual(values = c("0-5" = "#440154", "5-10" = "#3b528b",  
                               "10-20" = "#21918c", "20-30" = "#5ec962",  
                               "30-40" = "#fde725", "40+" = "#ffdf00"),  
                    name = "RRMSE") +  
  theme_minimal() +  
  labs(title = " ", fill = "Variable") +  
  theme(legend.position = c(0.9, 0.25)) +  
  xlim(-125, -65) + ylim(22.5, 50) +  
  annotation_scale(location = "bl", width_hint = 0.2)

