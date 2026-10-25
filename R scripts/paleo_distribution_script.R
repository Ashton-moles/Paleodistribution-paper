############################################################################################################
# Do oxygen and temperature control the distribution of marine ectotherms at the end Permian/early Triassic?
############################################################################################################

#Set Working directory
setwd("C:/Users/am14g22/OneDrive - University of Southampton/Paleodistribution paper")

#Load packages
install.packages("dismo")
install.packages("maps")
install.packages("raster")
install.packages("Rtools")
install.packages("installr")
install.packages("rJava")
install.packages("disdat")
install.packages("SDMtune")
install.packages("virtualspecies")
install.packages("ggpattern")
install.packages("ggridges")
install.packages("rayshader")
install.packages("forcats")
install.packages("magick")
install.packages("ggpubr")
install.packages("gridGraphics")
install.packages("ggtext")
install.packages("exactRankTests")
install.packages("PaleoClimR")
install.packages("pals")
install.packages("AquaEnv")
install.packages("paletteer")
install.packages("FSA")
install.packages("rnaturalearthdata")
install.packages("rgplates")
install.packages("palaeoverse")
install.packages("Cairo")
install.packages("PaleoClimR")
install.packages("RNetCDF")
install.packages("stars")

library(rayshader)
library(tidyr)
library(here)
library(broom)
library(FSA)
library(cowplot)
library(exactRankTests)
library(ggtext)
library(gridGraphics)
library(sf)
library(dismo)
library(ggpubr)
library(maps)
library(raster)
library(installr)
library(ggpattern)
library(rJava)
library(disdat)
library(SDMtune)
library(ggplot2)
library(paletteer)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgplates)
library(tidyverse)
library(palaeoverse)
library(Cairo)
library(PaleoClimR)
library(sp)
library(RNetCDF)
library(stars)
library(sp)
library(terra)
library(virtualspecies)
library(ggridges)
library(forcats)
library(patchwork)
library(magick)
library(corrplot)
library(AquaEnv)
library(pals)
devtools::install_github("richardstockey/PaleoClimR", force = TRUE)

# Load data
load("data/example_epme_cleaned_binned.RData")

# Create output folder
dir.create("Plots", showWarnings = FALSE)
dir.create("Statistics", showWarnings = FALSE)
dir.create("Data", showWarnings = FALSE)
dir.create("ENM_prediction_plots_pre", showWarnings = FALSE)
dir.create("response_curves_pre", showWarnings = FALSE)
dir.create("ENM_prediction_plots_post", showWarnings = FALSE)
dir.create("Response_curves_post", showWarnings = FALSE)


################################################################################
# Adding new columns and cleaning existing ones
################################################################################

# Add age column
stage_ages <- c(
  "48" = 264.28, 
  "49" = 259.51,   
  "50" = 254.14,  
  "51" = 251.90, 
  "52" = 251.2,  
  "53" = 247.2,  
  "54" = 242
  )

epme_data$age <- stage_ages[as.character(epme_data$stg)]

# Paleorotate data using PALEOMAP model
epme_data <- palaeorotate(
  occdf = epme_data,
  lng = "lng",
  lat = "lat",
  age = "age",
  model = "PALEOMAP",
  method = "point",
  uncertainty = FALSE,
  round = 1
)

# Delete old paleocords and rename columns
epme_data <- epme_data[,-46:-47]
epme_data <- rename(epme_data, paleolat = p_lat)
epme_data <- rename(epme_data, paleolng = p_lng)
epme_data <- rename(epme_data, species = accepted_name)

# Get rid of terrestrial taxa
epme_data <- epme_data %>% filter(!taxon_environment %in% c("freshwater","terrestrial", "brackish,freshwater,terrestrial", ""))

# Cleaning diets
epme_data <- epme_data %>% 
  mutate(diet = case_when(diet == "" ~ "unknown", is.na(diet) ~ "unknown",
 
diet %in% c("microcarnivore, carnivore",
             "durophage, carnivore",
             "microcarnivore",
             "carnivore, piscivore",
             "microcarnivore, photosymbiotic",
             "carnivore, detritivore") ~ "carnivore",
 
 diet %in% c("detritivore, suspension feeder",
             "suspension feeder, carnivore",
             "suspension feeder, photosymbiotic",
             "carnivore, suspension feeder",
             "suspension feeder, detritivore",
             "deposit feeder, suspension feeder",
             "photosymbiotic, suspension feeder") ~ "suspension feeder",
                                                   
 diet %in% c("detritivore, grazer",
             "grazer, deposit feeder",
             "grazer, detritivore",
             "grazer, photosymbiotic") ~ "grazer",
 
 diet %in% c("browser, omnivore",
             "omnivore, grazer",
             "grazer, omnivore") ~ "omnivore",
 
 diet %in% c("deposit feeder, chemosymbiotic") ~ "deposit feeder",
 
 diet %in% c("chemosymbiotic",
             "parasite",
             "durophage",
             "coprophage",
             "detritivore",
             "herbivore",
             "grazer",
             "deposit feeder") ~ "other",
 
 TRUE ~ diet))


# Cleaning vision
epme_data <- epme_data %>% mutate(vision = case_when(vision == "" ~ "unknown",
                                                     is.na(vision) ~ "unknown",
                                                     TRUE ~ vision))



# Cleaning life habits
epme_data <- epme_data %>% 
  mutate(life_habit = case_when(
    life_habit == "" ~ "unknown",
    is.na(life_habit) ~ "unknown",
    life_habit %in% c("low-level epifaunal", 
                      "intermediate-level epifaunal", 
                      "upper-level epifaunal", 
                      "epifaunal, gregarious",
                      "epifaunal, solitary, polymorph, depth=surface",
                      "intermediate-level epifaunal, colonial, clonal",
                      "epifaunal, colonial, clonal",
                      "epifaunal, gregarious, clonal") ~ "epifaunal",
    
    life_habit %in% c("deep infaunal", 
                      "semi-infaunal", 
                      "shallow infaunal") ~ "infaunal",
    
    life_habit == "nektobenthic, solitary" ~ "nektobenthic",
    
    life_habit %in% c("solitary", "amphibious", "aquatic", "boring", "colonial, clonal", "planktonic") ~ "other",
    
    TRUE ~ life_habit
    
  ))

# Clean motility
epme_data <- epme_data %>% mutate(motility = recode(motility,
                                                    "stationary, attached" = "stationary",
                                                    "stationary, attached, epibiont" = "stationary",
                                                    "fast-moving" = "actively mobile",
                                                    "slow-moving" = "actively mobile",
                                                    "facultatively mobile, attached" = "facultatively mobile",
                                                    "facultatively mobile, epibiont" = "epibiont",
                                                    "passively mobile, epibiont" = "epibiont"))

epme_data <- epme_data %>%
  filter(motility != "" & !is.na(motility))

write.csv(epme_data, here("Data", "epme_data.csv"))

################################################################################
# Seperating the data 
################################################################################

epme_data <- epme_data %>% 
  filter(!is.na(paleolat))

# Filter for distinct locations
epme_data <- epme_data %>% 
  distinct(genus, paleolng, paleolat, stg, .keep_all = TRUE)

# pre and post
epme_data_pre <- epme_data %>% filter(stg <= 51)
epme_data_post <- epme_data %>% filter(stg > 51)

# Separately list genera who were found before the extinction event and after
global_preextinction <- unique(epme_data$genus[epme_data$stg <= 51])
global_postextinction <- unique(epme_data$genus[epme_data$stg > 51])

# Find genera who survived and seperate to pre/post
global_survivors <- intersect(global_preextinction, global_postextinction)
global_survivor_epme_data <- epme_data %>% filter(genus %in% global_survivors)
global_survivor_epme_data_pre <- global_survivor_epme_data %>% filter(stg <= 51)
global_survivor_epme_data_post <- global_survivor_epme_data %>% filter(stg > 51)


# Find genera that went extinct
global_dead_vect <- setdiff(global_preextinction, global_survivors)
global_dead <- epme_data %>% filter(genus %in% global_dead_vect)

# Find genera that only occurred after the extinction event
global_originators_vect <- setdiff(global_postextinction, global_preextinction)
global_originators <- epme_data %>% filter(genus %in% global_originators_vect)

# Create seperate df for coliniarity tests with all variables
epme_data_pre_colintest <- epme_data_pre
epme_data_post_colintest_post <- epme_data_post


########################################################################################################
# Paleoclimate data
########################################################################################################

# Load climate model experiment 1
exp1 <- "HulseNatGeo2020/HulseNatGeo2020/01_endPermian_background_02.5CO2_1.00PO4_Tmp_98sink_sulfE6"
exp2 <- "HulseNatGeo2020/HulseNatGeo2020/02_StartWarming_07.5CO2_1.50PO4_Tmp_98sink_sulfE6"
exp3 <- "HulseNatGeo2020/HulseNatGeo2020/03_MainExtinction_20.0CO2_2.00PO4_Tmp_98sink_sulfE6"


# Load environmental data for experiment 1 (pre extinciton)
lvl2_O2_1 <- cGENIE.data(var = "ocn_O2", experiment = exp1, depth.level = 2)
lvl2_temp_1 <- cGENIE.data(var = "ocn_temp", experiment = exp1, depth.level = 2)
lvl2_H2S_1 <- cGENIE.data(var = "ocn_H2S", experiment = exp1, depth.level = 2)
lvl2_pH_1 <- cGENIE.data(var = "misc_pH", experiment = exp1, depth.level = 2)
lvl2_sal_1 <- cGENIE.data(var = "ocn_sal", experiment = exp1, depth.level = 2)
lvl2_Ca_1 <- cGENIE.data(var = "ocn_Ca", experiment = exp1, depth.level = 2)
lvl2_ALK_1 <- cGENIE.data(var = "ocn_ALK", experiment = exp1, depth.level = 2)
lvl2_DIC_1 <- cGENIE.data(var = "ocn_DIC", experiment = exp1, depth.level = 2)
lvl2_Fe_1 <- cGENIE.data(var = "ocn_Fe", experiment = exp1, depth.level = 2)
lvl2_Mg_1 <- cGENIE.data(var = "ocn_Mg", experiment = exp1, depth.level = 2)
lvl2_PO4_1 <- cGENIE.data(var = "ocn_PO4", experiment = exp1, depth.level = 2)
lvl2_SO4_1 <- cGENIE.data(var = "ocn_SO4", experiment = exp1, depth.level = 2)
lvl2_L_1 <- cGENIE.data(var = "ocn_L", experiment = exp1, depth.level = 2)
lvl2_FeL_1 <- cGENIE.data(var = "ocn_FeL", experiment = exp1, depth.level = 2)
lvl2_DOM_1 <- cGENIE.data(var = "ocn_DOM_C", experiment = exp1, depth.level = 2)
lvl2_POC_1 <- cGENIE.data(var = "bio_fparttot_POC", experiment = exp1, depth.level = 2)

# Load environmental data for experiment 3 (post extinction)
lvl2_O2_3 <- cGENIE.data(var = "ocn_O2", experiment = exp3, depth.level = 2)
lvl2_temp_3 <- cGENIE.data(var = "ocn_temp", experiment = exp3, depth.level = 2)
lvl2_H2S_3 <- cGENIE.data(var = "ocn_H2S", experiment = exp3, depth.level = 2)
lvl2_pH_3 <- cGENIE.data(var = "misc_pH", experiment = exp3, depth.level = 2)
lvl2_sal_3 <- cGENIE.data(var = "ocn_sal", experiment = exp3, depth.level = 2)
lvl2_Ca_3 <- cGENIE.data(var = "ocn_Ca", experiment = exp3, depth.level = 2)
lvl2_ALK_3 <- cGENIE.data(var = "ocn_ALK", experiment = exp3, depth.level = 2)
lvl2_DIC_3 <- cGENIE.data(var = "ocn_DIC", experiment = exp3, depth.level = 2)
lvl2_Fe_3 <- cGENIE.data(var = "ocn_Fe", experiment = exp3, depth.level = 2)
lvl2_Mg_3 <- cGENIE.data(var = "ocn_Mg", experiment = exp3, depth.level = 2)
lvl2_PO4_3 <- cGENIE.data(var = "ocn_PO4", experiment = exp3, depth.level = 2)
lvl2_SO4_3 <- cGENIE.data(var = "ocn_SO4", experiment = exp3, depth.level = 2)
lvl2_L_3 <- cGENIE.data(var = "ocn_L", experiment = exp3, depth.level = 2)
lvl2_FeL_3 <- cGENIE.data(var = "ocn_FeL", experiment = exp3, depth.level = 2)
lvl2_DOM_3 <- cGENIE.data(var = "ocn_DOM_C", experiment = exp3, depth.level = 2)
lvl2_POC_3 <- cGENIE.data(var = "bio_fparttot_POC", experiment = exp3, depth.level = 2)


# Attempt at metabolic index
#' 
#' press <- function(d, lat)
#' {
#'   press <- c()
#'   for (i in 1:length(d))
#'   {
#'     press <- c(press, gauge_p(d=d[[i]], lat=lat[[i]]))
#'   }
#'   return(press)
#' }
#' 
#' p_corr_pp <- function(pp, p, t, species="O2")
#' {
#'   nu <- switch(species,
#'                O2 =mean(31.9, 31.9, 32.2, 32.3, 32, 31.7),      # cm^3/mol (molar volume of O2 here)
#'                CO2=mean(35.3, 34.4),
#'                N2 =mean(33.2, 33.2, 33.5),
#'                Ar =mean(32.3, 32.1),
#'                He =mean(29.7, 29.7))
#'   R  <- 8.314472  # J/(mol*K) = Nm/(mol*K) = kg*m^2/s^2/(mol*K)
#'   P1 <- 0         # atm = 1.01325 bar = 1.01325*10^5 N/m^2
#'   T  <- t+273.15  # K
#'   p  <- p/1.01325 # atm = 1.01325 bar = 1.01325*10^5 N/m^2
#'   
#'   return(exp(-(nu*(P1-p)/(R*T*10)))*pp)
#' }
#' 
#' pO2_ <- function(O2, t, S, p)
#' {
#'   thet    <- theta(S=S, t=t, p=p)      # deg C
#'   O2_sat  <- O2sat_new(S=S, t=thet)    # umolO2/kg-soln
#'   xO2air  <- 0.20946                   # mole fraction: dimensionless
#'   
#'   pO2  <- (O2/O2_sat)*(xO2air*(1-sat_vap_press(S=S, t=thet))*1e6)   # at surface: total pressure = 1 atm 
#'   attr(pO2, "unit") <- "uatm"
#'   
#'   return(pO2)
#' }
#' 
#' theta <- function(S, t, p)
#' {
#'   PR <- 0                          # reference pressure
#'   H  <- PR - p                     # pressure differential
#'   
#'   tpar  <- c(0.5, 0.29289322,   1.707106781, 0)
#'   qpar1 <- c(1,   0.58578644,   3.414213562, 0)
#'   qpar2 <- c(0,   0.121320344, -4.121320344, 1)
#'   
#'   Q <- 0; XK <- 0 
#'   for (i in 1:4)
#'   {
#'     XK <- H*delta_adiabat(S, t, p)
#'     
#'     t  <- t + tpar[[i]]*(XK - Q)
#'     Q  <- qpar1[[i]]*XK + qpar2[[i]]*Q
#'     
#'     if (!(i==2)) {p  <- p + 0.5*H}
#'   }
#'   theta = t + (XK - 2*Q)/6
#'   
#'   attr(theta, "unit") <- "degC"
#'   return(theta)
#' }
#' 
#' delta_adiabat <- function(S, t, p)
#' {
#'   DS <- S - 35    
#'   p  <- p*10                            # convert bar to dbar
#'   
#'   delta_adiabat <- ((((-2.1687e-16*t + 1.8676e-14)*t - 4.6206e-13)*p + ((2.7759e-12*t - 1.1351e-10)*DS
#'                                                                         +((-5.4481e-14*t + 8.733e-12)*t - 6.7795e-10)*t + 1.8741e-8))*p
#'                     +(-4.2393e-8*t + 1.8932e-6)*DS+((6.6228e-10*t - 6.836e-8)*t + 8.5258e-6)*t + 3.5803e-5)
#'   
#'   delta_adiabat <- delta_adiabat*10     # convert /dbar to /bar
#'   
#'   attr(delta_adiabat, "unit") <- "degC/bar"
#'   return(delta_adiabat)
#' }
#' 
#' O2sat_new <- function(S, t)
#' {
#'   #'CALCULATE OXYGEN CONCENTRATION AT SATURATION
#'   #'Source:  Oxygen solubility in seawater: Better fitting equations.
#'   #'          Garcia & Gordon (1992) L&O 37: 1307-1312.
#'   #'Input:       S = Salinity (pss-78)
#'   #'             T = Temp (deg C)
#'   #'Output:      Oxygen saturation at one atmosphere (umol/kg).
#'   
#'   #'DEFINE CONSTANTS, ETC FOR SATURATION CALCULATION
#'   
#'   #' The constants used are for units of umol O2/kg.
#'   T <- t
#'   
#'   A0 = 5.80871
#'   A1 = 3.20291
#'   A2 = 4.17887
#'   A3 = 5.10006
#'   A4 = -0.0986643
#'   A5 = 3.80369
#'   
#'   B0 = -0.00701577
#'   B1 = -0.00770028
#'   B2 = -0.0113864
#'   B3 = -0.00951519
#'   
#'   C0 = -0.000000275915
#'   
#'   #'CALCULATE Ts FROM T (deg C):
#'   
#'   TS = log((298.15 - T) / (273.15 + T))
#'   
#'   #'CALCULATE O2 SATURATION in umol/kg:
#'   
#'   A = ((((A5 * TS + A4) * TS + A3) * TS + A2) * TS + A1) * TS + A0
#'   B = ((B3 * TS + B2) * TS + B1) * TS + B0
#'   O2sat_new = exp(A + S * (B + S * C0))
#' }
#' 
#' sat_vap_press <- function(S, t)
#' {
#'   T = t + 273.15
#'   return(exp(24.4543 - (6745.09/T) - 4.8489*log(T/100) - 0.000544*S))
#' }
#' 
#' # colnames(lvl2_O2_1)[colnames(lvl2_O2_1) == 'var'] <- "O2"
#' # colnames(lvl2_sal_1)[colnames(lvl2_sal_1) == 'var'] <- "salinity"
#' # colnames(lvl2_temp_1)[colnames(lvl2_temp_1) == 'var'] <- "temperature" 
#' 
#' pO2_lvl1 <- left_join(lvl2_O2_1, lvl2_sal_1)
#' pO2_lvl1 <- left_join(pO2_lvl1, lvl2_temp_1)
#' pO2_lvl1$depth <- 200
#' 
#' pO2 <- function(O2, S, t, d, lat) 
#'   
#' { 
#'   
#'   p   <- press(d=d, lat=lat) 
#'   
#'   pO2 <- (p_corr_pp(pp=pO2_(O2=O2, t=t, S=S, p=p), p=p, t=t))/1000 
#'   
#'   attr(pO2, "unit") <- "matm" 
#'   
#'   return(pO2) 
#'   
#' }
#' 
#' pO2_lvl1$pO2 <- mapply(
#'   FUN  = pO2,
#'   O2   = pO2_lvl1$O2,
#'   S    = pO2_lvl1$sal,
#'   t    = pO2_lvl1$temp,
#'   d    = pO2_lvl1$depth,
#'   lat  = pO2_lvl1$lat.mid
#' )

# Plot occurrences, sanity check #
# Use var with no data to act as blank map 
lvl2_blank_1 <- cGENIE.data(var = "bio_fparttot_CaCO3", experiment = exp1, depth.level = 2)
lvl2_blank_map <- cGENIE.map(var = "bio_fparttot_CaCO3", experiment = exp1, depth.level = 2, dims = 3, min.value = 6, max.value = 8, intervals = .2, unit.factor = 1, scale.label = NA)

# Make as a sf object
sf <- st_as_sf(epme_data_pre, coords = c("paleolng", "paleolat"), crs = 4326)

# Plot
lvl2_blank_map <- lvl2_blank_map + geom_sf(data = sf, col = "black", size = 4, alpha = 1, pch = 20) +
  labs(title = "pre-extinction occurences")

ggsave(
  filename = "Plots/pre-extinction occurences plot.png",
  plot = lvl2_blank_map,
  width = 8,
  height = 6,
  dpi = 300
)

# Define datasets with all variables (pre and post)
            colintest <- list(
              pH = lvl2_pH_1, O2 = lvl2_O2_1, temp = lvl2_temp_1, H2S = lvl2_H2S_1,
              SO4 = lvl2_SO4_1,sal = lvl2_sal_1, ca = lvl2_Ca_1, POC = lvl2_POC_1,
              ALK = lvl2_ALK_1, DIC = lvl2_DIC_1, Fe = lvl2_Fe_1, Mg = lvl2_Mg_1,
              L = lvl2_L_1, FeL = lvl2_FeL_1
              )
            
            colintest_post <- list(
              pH = lvl2_pH_3, O2 = lvl2_O2_3, temp = lvl2_temp_3, H2S = lvl2_H2S_3, 
              POC = lvl2_POC_3, SO4 = lvl2_SO4_3, sal = lvl2_sal_3, ca = lvl2_Ca_3,
              ALK = lvl2_ALK_3, DIC = lvl2_DIC_3, Fe = lvl2_Fe_3, Mg = lvl2_Mg_3, 
              L = lvl2_L_3, FeL = lvl2_FeL_3
              )

# Define datasets with important variables (pre and post)
env_datasets <- list(
  pH = lvl2_pH_1, O2 = lvl2_O2_1, temp = lvl2_temp_1, H2S = lvl2_H2S_1, POC = lvl2_POC_1, SO4 = lvl2_SO4_1
)

env_datasets3 <- list(
  pH = lvl2_pH_3, O2 = lvl2_O2_3, temp = lvl2_temp_3, H2S = lvl2_H2S_3, POC = lvl2_POC_3, SO4 = lvl2_SO4_3
)

# Datasets without collniearity (pre - no o2 and no temp)
env_datasets_sans_o2 <- list(
  pH = lvl2_pH_1,temp = lvl2_temp_1, H2S = lvl2_H2S_1, POC = lvl2_POC_1, SO4 = lvl2_SO4_1)
env_datasets_sans_temp <- list(
  pH = lvl2_pH_1, O2 = lvl2_O2_1, H2S = lvl2_H2S_1, POC = lvl2_POC_1, SO4 = lvl2_SO4_1)

# Post no pH and no o2
env_datasets_sans_pH_post <- list(
  temp = lvl2_temp_3, O2 = lvl2_O2_3, H2S = lvl2_H2S_3, POC = lvl2_POC_3, SO4 = lvl2_SO4_3)
env_datasets_sans_o2_post <- list(
  pH = lvl2_pH_3, temp = lvl2_temp_3, H2S = lvl2_H2S_3, POC = lvl2_POC_3, SO4 = lvl2_SO4_3)

#################################################################################################################################
# Extract environmental data at occurrences
#################################################################################################################################

# Function to extract environmental data of the nearest ocean point if NA
nearest_ocean_extractor <- function(var_df, max_dist = 5) { 
  function(fossil_lon, fossil_lat) {
    
    if (is.na(fossil_lon) || is.na(fossil_lat)) {
      return(NA_real_)
    }
    distances <- sqrt(
      (var_df$lon.mid - fossil_lon)^2 +
        (var_df$lat.mid - fossil_lat)^2
    )
    valid_cells <- !is.na(var_df$var)
    
    if (!any(valid_cells)) {
      return(NA_real_)
    }
    min_dist <- min(distances[valid_cells])
    
    if (min_dist > max_dist) {
      return(NA_real_)
    }
    valid_values <- var_df$var[valid_cells]
    valid_values[which.min(distances[valid_cells])]
  }
}
# Rasterise
template_rast <- rast(
  xmin = -180, xmax = 180,
  ymin = -90, ymax = 90,
  resolution = c(10, 10),
  crs = "EPSG:4326"
)

cGENIE_df_to_raster  <- function(
    df,
    template_rast,
    value_col = "var",
    lon_min = "lon.min",
    lon_max = "lon.max", 
    lat_min = "lat.min",
    lat_max = "lat.max",
    drop_na = FALSE,
    raster_fun = "mean"
) {
  
  stopifnot(
    all(c(value_col, lon_min, lon_max, lat_min, lat_max) %in% names(df))
  )
  
  if (drop_na) {
    df <- df |> filter(!is.na(.data[[value_col]]))
  }
  
  sf_cells <- df |>
    rowwise() |>
    mutate(
      geometry = list(st_polygon(list(matrix(
        c(.data[[lon_min]], .data[[lat_min]],
          .data[[lon_max]], .data[[lat_min]],
          .data[[lon_max]], .data[[lat_max]],
          .data[[lon_min]], .data[[lat_max]],
          .data[[lon_min]], .data[[lat_min]]),
        ncol = 2,
        byrow = TRUE
      ))))
    ) |>
    ungroup() |>
    st_as_sf(crs = 4326)
  
  r <- rasterize(
    vect(sf_cells),
    template_rast,
    field = value_col,
    fun = raster_fun,
    touches = FALSE
  )
}

tolerate_dist <- 5

###################################
# Extracting variable data (*pre*)
###################################
for (var_name in names(env_datasets)) {
  message(paste("Extracting:", var_name))
  
  extractor <- nearest_ocean_extractor(env_datasets[[var_name]], max_dist = tolerate_dist)
  epme_data_pre[[var_name]] <- sapply(1:nrow(epme_data_pre), function(i) {
    extractor(epme_data_pre$paleolng[i], epme_data_pre$paleolat[i])
  })
}

# Remove NAs from points 
epme_data_pre <- epme_data_pre %>% 
  filter(!is.na(temp))

# Generate rasters
raster_list <- lapply(env_datasets, function(df) {
  cGENIE_df_to_raster(
    df = df,
    value_col = "var",
    template_rast = template_rast
  )
})

# Convert the list to a multi-layer SpatRaster object
Env_rasters_all_pre <- rast(raster_list)
Env_rasters_all_pre <- mask(Env_rasters_all_pre, sum(Env_rasters_all_pre))

###################################
# Extracting variable data (*post*)
###################################
for (var_name in names(env_datasets3)) {
  message(paste("Extracting:", var_name))
  
  extractor <- nearest_ocean_extractor(env_datasets3[[var_name]], max_dist = tolerate_dist)
  epme_data_post[[var_name]] <- sapply(1:nrow(epme_data_post), function(i) {
    extractor(epme_data_post$paleolng[i], epme_data_post$paleolat[i])
  })
}

# Remove NAs from points 
epme_data_post <- epme_data_post %>% 
  filter(!is.na(temp))

# Generate rasters
raster_list <- lapply(env_datasets3, function(df) {
  cGENIE_df_to_raster(
    df = df,
    value_col = "var",
    template_rast = template_rast
  )
})

# Convert the list to a multi-layer SpatRaster object
Env_rasters_all_post <- rast(raster_list)
Env_rasters_all_post <- mask(Env_rasters_all_post, sum(Env_rasters_all_post))

#################################################
# Extracting variable data for colin test (*pre*)
#################################################
for (var_name in names(colintest)) {
  message(paste("Extracting:", var_name))
  
  extractor <- nearest_ocean_extractor(colintest[[var_name]], max_dist = tolerate_dist)
  epme_data_pre_colintest[[var_name]] <- sapply(1:nrow(epme_data_pre_colintest), function(i) {
    extractor(epme_data_pre_colintest$paleolng[i], epme_data_pre_colintest$paleolat[i])
  })
}

# Remove NAs from points 
epme_data_pre_colintest <- epme_data_pre_colintest %>% 
  filter(!is.na(temp))

# Generate rasters
raster_list <- lapply(colintest, function(df) {
  cGENIE_df_to_raster(
    df = df,
    value_col = "var",
    template_rast = template_rast
  )
})

# Convert the list to a multi-layer SpatRaster object
Env_rasters_colintest <- rast(raster_list)
Env_rasters_colintest <- mask(Env_rasters_colintest, sum(Env_rasters_colintest))

#################################################
# Extracting variable data for colin test (*post*)
#################################################
for (var_name in names(colintest_post)) {
  message(paste("Extracting:", var_name))
  
  extractor <- nearest_ocean_extractor(colintest_post[[var_name]], max_dist = tolerate_dist)
  epme_data_post_colintest_post[[var_name]] <- sapply(1:nrow(epme_data_post_colintest_post), function(i) {
    extractor(epme_data_post_colintest_post$paleolng[i], epme_data_post_colintest_post$paleolat[i])
  })
}

# Remove NAs from points 
epme_data_post_colintest_post <- epme_data_post_colintest_post %>% 
  filter(!is.na(temp))

# Generate rasters
raster_list <- lapply(colintest_post, function(df) {
  cGENIE_df_to_raster(
    df = df,
    value_col = "var",
    template_rast = template_rast
  )
})

# Convert the list to a multi-layer SpatRaster object
Env_rasters_colintest_post <- rast(raster_list)
Env_rasters_colintest_post <- mask(Env_rasters_colintest_post, sum(Env_rasters_colintest_post))

##########################################
# Run collinearity tests
##########################################

# (*pre*)
raster_legacy_list_colintest <- lapply(as.list(Env_rasters_colintest), raster)
worldclim_stack_colintest <- stack(raster_legacy_list_colintest)

png(filename = "Statistics/Pre-extinction Collinearity dendrogram.png", width = 8, height = 6, units = "in", res = 300)
cor_matrix_pre <- removeCollinearity(worldclim_stack_colintest, plot = TRUE, multicollinearity.cutoff = 0.7, method = "pearson")
dev.off()

stack_cor <- layerStats(worldclim_stack_colintest, 'pearson', na.rm=TRUE)
cor_mat <- stack_cor$`pearson correlation coefficient`
p_mat <- cor_mat
cross_out_mat <- ifelse(abs(cor_mat) > 0.7, 1, 0)
diag(cross_out_mat) <- 0
no_cross_lower <- matrix(0, nrow = nrow(cor_mat), ncol = ncol(cor_mat))

# (*post*)
raster_legacy_list_colintest_post <- lapply(as.list(Env_rasters_colintest_post), raster)
worldclim_stack_colintest_post <- stack(raster_legacy_list_colintest_post)

png(filename = "Statistics/Post-extinction Collinearity dendrogram.png", width = 8, height = 6, units = "in", res = 300)
cor_matrix_post <- removeCollinearity(worldclim_stack_colintest_post, plot = TRUE, multicollinearity.cutoff = 0.7, method = "pearson")
dev.off()

stack_cor_post <- layerStats(worldclim_stack_colintest_post, 'pearson', na.rm=TRUE)
cor_matrix_post <- stack_cor_post$`pearson correlation coefficient`

cross_out_matrix_post <- ifelse(abs(cor_matrix_post) > 0.7, 1, 0)

png(filename = "statistics/Pre and post extinction collinearity heatmaps.png",  width = 8, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

cor_plot_pre <- corrplot.mixed(cor_mat, 
                               upper = "circle",
                               lower = "number",
                               p.mat = cross_out_mat, 
                               sig.level = 0.5,
                               insig = "pch",       
                               pch = "X",
                               pch.col = "black",
                               pch.cex = 1.5,
                               lower.p.mat = no_cross_lower, 
                               tl.pos = "lt",      
                               tl.col = "black",
                               lower.col = "black",
                               number.cex = 0.7,
                               cl.lim = c(-1, 1),
                               mar = c(0, 0, 2, 0)) 

mtext("Pre-Extinction", side = 3, line = 1, adj = 0.5, cex = 1.2, font = 2) 
  
  corrplot(cor_mat, 
          add = TRUE,             
          type = "lower", 
          method = "number", 
          col = "black",          
          diag = FALSE, 
          tl.pos = "n",           
          cl.pos = "n",           
          number.cex = 0.8)

# (*post*)

cor_plot_post <- corrplot.mixed(cor_matrix_post, 
               upper = "circle",
               lower = "number",
               p.mat = cross_out_matrix_post, 
               sig.level = 0.5,
               insig = "pch",       
               pch = "X",
               pch.col = "black",
               pch.cex = 1.5,
               tl.pos = "lt",      
               tl.col = "black",
               lower.col = "black",
               number.cex = 0.7,
               cl.lim = c(-1, 1),mar = c(0, 0, 2, 0)) 

  mtext("Post-Extinction", side = 3, line = 1, adj = 0.5, cex = 1.2, font = 2) 

corrplot(cor_matrix_post, 
         add = TRUE,             
         type = "lower", 
         method = "number", 
         col = "black",          
         diag = FALSE, 
         tl.pos = "n",           
         cl.pos = "n",           
         number.cex = 0.8)

cor_plot_combined <- recordPlot()
dev.off()

#################################################################################################
# Plot EMNs and save maxent results
#################################################################################################

# Run maxent (*pre*)
genera_pre <- unique(epme_data_pre$genus)
results_list_pre <- list()

# Convert the terra SpatRaster into a list of old 'Raster' objects
raster_legacy_list2 <- lapply(as.list(Env_rasters_all_pre), raster)

# Stack them into a RasterStack
worldclim_stack <- stack(raster_legacy_list2)
env_vars <- names(worldclim_stack)
names(worldclim_stack) <- names(Env_rasters_all_pre)

cor_matrix <- removeCollinearity(worldclim_stack, plot = TRUE, multicollinearity.cutoff = 0.7)

# Loop over all genera 
for (i in 1:length(genera_pre)) {  
  genus_name <- genera_pre[i]
  var_data <- epme_data_pre %>% filter(genus == genus_name)
  
  if (nrow(var_data) < 10) next
  
  varpoints <- as.matrix(var_data[, c('paleolng', 'paleolat')])
  group <- kfold(var_data, 10)
  pres_train <- varpoints[group != 1, ]
  pres_test <- varpoints[group == 1, ]  
  
  tryCatch({
    # Generate background points
    backg2 <- randomPoints(worldclim_stack, n=10000)
    colnames(backg2) <- c("lon","lat")
    group_backg2 <- kfold(backg2, 10)  
    backg_train2 <- backg2[group_backg2 != 1, ]
    backg_test2 <- backg2[group_backg2 == 1, ]
    
    # Run maxent
    xm <- maxent(worldclim_stack, pres_train, args = c("noproduct", "nothreshold", "nohinge", "betamultiplier=2.0", "jackknife=TRUE"))
    
    png(filename = paste0("response_curves_pre/", genus_name, "_response_curves.png"),
        width = 10, height = 8, units = "in", res = 300)
    response(xm)
    dev.off()
    
    eval_result <- evaluate(pres_test, backg_test2, xm, worldclim_stack)

    # Tss
    thr_table <- threshold(eval_result)
    best_thr <- thr_table$spec_sens
    idx <- which(eval_result@t == best_thr)
    current_tss <- eval_result@TPR[idx] + eval_result@TNR[idx] - 1
    
    #Plot prediction
    p <- predict(worldclim_stack, xm, progress='')

    png(filename = paste0("ENM_prediction_plots_pre/", genus_name, "_distribution_pre.png"),
        width = 12, height = 5, units = "in", res = 300)

    plot(p, main = paste(genus_name, "Pre-Extinction Habitat- Pre-Extinction\nAUC:",
                         round(eval_result@auc, 3)))
    points(pres_train, pch = 20, cex = 0.5, col = "black")

    dev.off()
    
    # Extract all maxent results
    all_results_pre <- xm@results[, 1]
    names(all_results_pre) <- rownames(xm@results)
    
    # Extract evaluation metrics
    eval_metrics <- c(
      test_AUC = eval_result@auc,
      test_TSS = current_tss,
      test_cor = eval_result@cor,
      test_pcor = eval_result@pcor,
      np = eval_result@np,  
      na = eval_result@na   
    )
    
    # Combine everything into one data frame row
    results_list_pre[[genus_name]] <- data.frame(
      genus = genus_name,
      t(all_results_pre),
      t(eval_metrics)
    )
  }, error = function(e) NULL)
}

# Save as data frame and csv
ENM_results_pre <- bind_rows(results_list_pre)
write.csv(ENM_results_pre,"Data/ENM_results_pre_all.csv", row.names = FALSE)

#####################
# Run maxent (*post*)
#####################

genera_post <- unique(epme_data_post$genus)
results_list_post <- list()

# Convert the terra SpatRaster into a list of old 'Raster' objects
raster_legacy_list3 <- lapply(as.list(Env_rasters_all_post), raster)

# Stack them into a RasterStack
worldclim_stack3 <- stack(raster_legacy_list3)
names(worldclim_stack3) <- names(Env_rasters_all_post)

cor_matrix <- removeCollinearity(worldclim_stack3, plot = TRUE, multicollinearity.cutoff = 0.7) 


for (i in 1:length(genera_post)) { 
  genus_name <- genera_post[i]
  var_data <- epme_data_post %>% filter(genus == genus_name)
  
  if (nrow(var_data) < 10) next
  
  varpoints <- as.matrix(var_data[, c('paleolng', 'paleolat')])
  group <- kfold(var_data, 10)
  pres_train <- varpoints[group != 1, ]
  pres_test <- varpoints[group == 1, ]  
  
  tryCatch({
    # Generate background points
    backg3 <- randomPoints(worldclim_stack3, n=10000)
    colnames(backg3) <- c("lon","lat")
    group_backg3 <- kfold(backg3, 10)  
    backg_train3 <- backg3[group_backg3 != 1, ]
    backg_test3 <- backg3[group_backg3 == 1, ]
    
    # Run maxent
    xm <- maxent(worldclim_stack3, pres_train, args = c("noproduct", "nothreshold", "nohinge", "betamultiplier=2.0", "jackknife=TRUE"))
    
    png(filename = paste0("response_curves_post/", genus_name, "_response_curves.png"),
        width = 10, height = 8, units = "in", res = 300)
    response(xm)
    dev.off()
    
    
    # Evaluate model
    eval_result <- evaluate(pres_test, backg_test3, xm, worldclim_stack3)

    # Tss
    thr_table <- threshold(eval_result)
    best_thr <- thr_table$spec_sens
    idx <- which(eval_result@t == best_thr)
    current_tss <- eval_result@TPR[idx] + eval_result@TNR[idx] - 1
    
    # Plot prediction
    p <- predict(worldclim_stack3, xm, progress='')
    png(filename = paste0("ENM_prediction_plots_post/", genus_name, "_distribution_post.png"),
        width = 12, height = 5, units = "in", res = 300)

    plot(p, main = paste(genus_name, "Post-Extinction Habitat- Post-Extinction\nAUC:",
                         round(eval_result@auc, 3)))
    points(pres_train, pch = 20, cex = 0.5, col = "black")

    dev.off()

    # Extract all maxent results
    all_results_post <- xm@results[, 1]
    names(all_results_post) <- rownames(xm@results)
    
    # Extract evaluation metrics
    eval_metrics <- c(
      test_AUC = eval_result@auc,
      test_TSS = current_tss,
      test_cor = eval_result@cor,
      test_pcor = eval_result@pcor,
      np = eval_result@np,
      na = eval_result@na
    )
    
    # Combine everything into one data frame row
    results_list_post[[genus_name]] <- data.frame(
      genus = genus_name,
      t(all_results_post),
      t(eval_metrics)
      )
    
  }, error = function(e) NULL)
}

# Save as data frame and csv
ENM_results_post <- bind_rows(results_list_post)
write.csv(ENM_results_post,"Data/ENM_results_post.csv", row.names = FALSE)

########################
# Add columns to results
########################

# Add catagories to enm results and filter out NAs
genus_and_life_history_pre <- epme_data_pre %>% select(21:25, 119:122)
genus_and_life_history_pre <- unique(genus_and_life_history_pre)
genus_and_life_history_post<- epme_data_post %>% select(21:25, 119:122)
genus_and_life_history_post <- unique(genus_and_life_history_post)


ENM_results_pre <- ENM_results_pre %>% left_join(genus_and_life_history_pre, by = "genus")
ENM_results_post <- ENM_results_post %>% left_join(genus_and_life_history_post, by = "genus")

ENM_results_pre <- ENM_results_pre %>% filter(!is.na(test_cor.cor))
ENM_results_post <- ENM_results_post %>% filter(!is.na(test_cor.cor))

############################################################################################################
# EMN with collinear vars removed
############################################################################################################

process_env_set <- function(env_list, occ_data, template) {
  # Extract values to occurrences
  for (var_name in names(env_list)) {
    extractor <- nearest_ocean_extractor(env_list[[var_name]], max_dist = tolerate_dist)
    occ_data[[var_name]] <- sapply(1:nrow(occ_data), function(i) {
      extractor(occ_data$paleolng[i], occ_data$paleolat[i])
    })
  }
  
  # Filter NAs 
  first_var <- names(env_list)[1]
  occ_data <- occ_data %>% filter(!is.na(!!sym(first_var)))
  
  # Create Raster Stack
  r_list <- lapply(env_list, function(df) {
    cGENIE_df_to_raster(df = df, value_col = "var", template_rast = template)
  })
  r_stack <- rast(r_list)
  r_stack <- mask(r_stack, sum(r_stack))
  
  return(list(occ = occ_data, stack = r_stack))
}


pre_sans_o2_results <- process_env_set(env_datasets_sans_o2, epme_data_pre, template_rast)
epme_data_pre_sans_o2 <- pre_sans_o2_results$occ
Env_rasters_pre_sans_o2 <- pre_sans_o2_results$stack

pre_sans_temp_results <- process_env_set(env_datasets_sans_temp, epme_data_pre, template_rast)
epme_data_pre_sans_temp <- pre_sans_temp_results$occ
Env_rasters_pre_sans_temp <- pre_sans_temp_results$stack

post_sans_o2_results <- process_env_set(env_datasets_sans_o2_post, epme_data_post, template_rast)
epme_data_post_sans_o2 <- post_sans_o2_results$occ
Env_rasters_post_sans_o2 <- post_sans_o2_results$stack

post_sans_pH_results <- process_env_set(env_datasets_sans_pH_post, epme_data_post, template_rast)
epme_data_post_sans_pH <- post_sans_pH_results$occ
Env_rasters_post_sans_pH <- post_sans_pH_results$stack


# Run the ENM
scenarios <- list(
  "pre_sans_o2"   = list(stack = Env_rasters_pre_sans_o2,   occ = epme_data_pre_sans_o2),
  "pre_sans_temp" = list(stack = Env_rasters_pre_sans_temp, occ = epme_data_pre_sans_temp),
  "post_sans_o2"  = list(stack = Env_rasters_post_sans_o2, occ = epme_data_post),
  "post_sans_ph"  = list(stack = Env_rasters_post_sans_pH, occ = epme_data_post)
)

for (scen_name in names(scenarios)) {
  
  message(paste("--- Running Scenario:", scen_name, "---"))
  
  current_stack_raw <- scenarios[[scen_name]]$stack
  current_occ_data  <- scenarios[[scen_name]]$occ
  
  raster_legacy_list <- lapply(as.list(current_stack_raw), raster)
  worldclim_stack    <- stack(raster_legacy_list)
  
  genera_list <- unique(current_occ_data$genus)
  results_list <- list()
  
  for (i in 1:length(genera_list)) {
    genus_name <- genera_list[i]
    var_data <- current_occ_data %>% filter(genus == genus_name)
    
    if (nrow(var_data) < 10) next
    
    varpoints <- as.matrix(var_data[, c('paleolng', 'paleolat')])
    group <- kfold(var_data, 10)
    pres_train <- varpoints[group != 1, ]
    pres_test  <- varpoints[group == 1, ]  
    
    tryCatch({
      backg_points <- randomPoints(worldclim_stack, n=10000)
      colnames(backg_points) <- c("lon","lat")
      group_backg <- kfold(backg_points, 10)  
      backg_test  <- backg_points[group_backg == 1, ]
      
      xm <- maxent(worldclim_stack, pres_train, 
                   args = c("noproduct", "nothreshold", "nohinge", "betamultiplier=2.0", "jackknife=TRUE"))
      
      eval_result <- evaluate(pres_test, backg_test, xm, worldclim_stack)
      
      thr_table <- threshold(eval_result)
      best_thr  <- thr_table$spec_sens 
      idx       <- which(eval_result@t == best_thr)
      current_tss <- eval_result@TPR[idx] + eval_result@TNR[idx] - 1
      
      all_results <- xm@results[, 1]
      
      results_list[[genus_name]] <- data.frame(
        genus = genus_name,
        scenario = scen_name,
        t(all_results),
        test_AUC = eval_result@auc,
        test_TSS = current_tss,
        test_cor = eval_result@cor
      )
      
    }, error = function(e) NULL)
  }
  
  final_df_name <- paste0("ENM_results_", scen_name)
  assign(final_df_name, bind_rows(results_list), envir = .GlobalEnv)
  
  write.csv(get(final_df_name), paste0("Data/",final_df_name, ".csv"), row.names = FALSE)
}


##################################################################################################################
# Generate plots - context plots
#################################################################################################################

#################
# Occurrence maps
#################

# Use var with no data to act as blank map 
lvl2_blank_1 <- cGENIE.data(var = "bio_fparttot_CaCO3", experiment = exp1, depth.level = 2)
lvl2_blank_map <- cGENIE.map(var = "bio_fparttot_CaCO3", experiment = exp1, depth.level = 2, dims = 3, min.value = 6, max.value = 8, intervals = .2, unit.factor = 1, scale.label = NA)

# Make as a sf object
sf_pre <- st_as_sf(epme_data_pre, coords = c("paleolng", "paleolat"), crs = 4326)
sf_post <- st_as_sf(epme_data_post, coords = c("paleolng", "paleolat"), crs = 4326)

# Plot
lvl2_blank_map_pre <- lvl2_blank_map + geom_sf(data = sf_pre, col = "red", size = 4, alpha = 1, pch = 20)
lvl2_blank_map_post <- lvl2_blank_map + geom_sf(data = sf_post, col = "red", size = 5, alpha = 0.3, pch = 20)

# Stitch together
occurrence_maps <- (lvl2_blank_map_pre | lvl2_blank_map_post) + 
                    plot_annotation(
                      title = "occurrencees of end-Permian ectotherms",
                      subtitle = "Pre-extinction(left) vs Post-extinction(right)"
                      )
ggsave("Plots/occurrence_maps.png", width = 15, height = 19, dpi = 600)



####################
# Climate model maps
####################

lvl2_O2_map <- cGENIE.map(var = "ocn_O2", experiment = exp1, depth.level = 2, dims = 3, min.value = 0, max.value = 300, intervals = 30, unit.factor = 1e6, scale.label = "Oxygen (\u00B5mol/kg)", palette_name = paletteer::paletteer_c("grDevices::PuBuGn", 4, direction = 1)) + labs(title = "Pre-extinction")
lvl2_temp_map <- cGENIE.map(var = "ocn_temp", experiment = exp1, depth.level = 2, dims = 3, min.value = 0, max.value = 32, intervals = 2, unit.factor = 1, scale.label = "Temp (°C)", palette_name = paletteer::paletteer_c("viridis::inferno", 30, direction = 1))
lvl2_pH_map <- cGENIE.map(var = "misc_pH", experiment = exp1, depth.level = 2, dims = 3, min.value = 6.4, max.value = 7.6, intervals = 0.1, unit.factor = 1, scale.label = "pH", palette_name = paletteer::paletteer_c("grDevices::Reds", 30, direction = 1))
lvl2_DOM_map <- cGENIE.map(var = "ocn_DOM_C", experiment = exp1, depth.level = 2, dims = 3, min.value = 0, max.value = 45, intervals = 5, unit.factor = 1e6, scale.label = "DOM: carbon (mol kg-1)", palette_name = paletteer::paletteer_c("grDevices::Greens", 30, direction = -1))
lvl2_H2S_map <- cGENIE.map(var = "ocn_H2S", experiment = exp1, depth.level = 2, dims = 3, min.value = 0, max.value = 500, intervals = 50, unit.factor = 1e6, scale.label = "H2S (mol kg-1)", palette_name = paletteer::paletteer_d("ggsci::amber_material", 10, direction = 1))
lvl2_SO4_map <- cGENIE.map(var = "ocn_SO4", experiment = exp1, depth.level = 2, dims = 3, min.value = 0.027, max.value = 0.035, intervals = 0.002, unit.factor = 1, scale.label = "SO4 (mol kg-1)", palette_name = paletteer::paletteer_c("grDevices::RdPu", 30, direction = -1))

final_comparison <- lvl2_O2_map + lvl2_temp_map + lvl2_pH_map + lvl2_DOM_map + lvl2_H2S_map + lvl2_SO4_map

terra::writeRaster(exp1$ocn_O2, "Oxygen_Map.tif")

# post
lvl2_O2_map_3 <- cGENIE.map(var = "ocn_O2", experiment = exp3, depth.level = 2, dims = 3, min.value = 0, max.value = 300, intervals = 30, unit.factor = 1e6, scale.label = "Oxygen (\u00B5mol/kg)", palette_name = paletteer::paletteer_c("grDevices::PuBuGn", 4, direction = 1)) + labs(title = "Post-extinction")
lvl2_temp_map_3 <- cGENIE.map(var = "ocn_temp", experiment = exp3, depth.level = 2, dims = 3, min.value = 0, max.value = 32, intervals = 2, unit.factor = 1, scale.label = "Temp (°C)", palette_name = paletteer::paletteer_c("viridis::inferno", 30, direction = 1))
lvl2_pH_map_3 <- cGENIE.map(var = "misc_pH", experiment = exp3, depth.level = 2, dims = 3, min.value = 6.4, max.value = 7.6, intervals = 0.1, unit.factor = 1, scale.label = "pH", palette_name = paletteer::paletteer_c("grDevices::Reds", 30, direction = 1))
lvl2_DOM_map_3 <- cGENIE.map(var = "ocn_DOM_C", experiment = exp3, depth.level = 2, dims = 3, min.value = 0, max.value = 45, intervals = 5, unit.factor = 1e6, scale.label = "DOM: carbon (mol kg-1)", palette_name = paletteer::paletteer_c("grDevices::Greens", 30, direction = -1))
lvl2_H2S_map_3 <- cGENIE.map(var = "ocn_H2S", experiment = exp3, depth.level = 2, dims = 3, min.value = 0, max.value = 500, intervals = 50, unit.factor = 1e6, scale.label = "H2S (mol kg-1)", palette_name = paletteer::paletteer_d("ggsci::amber_material", 10, direction = 1))
lvl2_SO4_map_3 <- cGENIE.map(var = "ocn_SO4", experiment = exp3, depth.level = 2, dims = 3, min.value = 0.027, max.value = 0.035, intervals = 0.002, unit.factor = 1, scale.label = "SO4 (mol kg-1)", palette_name = paletteer::paletteer_c("grDevices::RdPu", 30, direction = -1))

final_comparison_3 <- lvl2_O2_map_3 + lvl2_temp_map_3 + lvl2_pH_map_3 + lvl2_DOM_map_3 + lvl2_H2S_map_3 + lvl2_SO4_map_3

# Stitch together
o2_maps <- lvl2_O2_map + lvl2_O2_map_3
temp_maps <- lvl2_temp_map + lvl2_temp_map_3
pH_maps <- lvl2_pH_map + lvl2_pH_map_3
DOM_maps <- lvl2_DOM_map + lvl2_DOM_map_3
H2S_maps <- lvl2_H2S_map + lvl2_H2S_map_3
SO4_maps <- lvl2_SO4_map + lvl2_SO4_map_3

final_comparison_both <- o2_maps / temp_maps / pH_maps / DOM_maps / H2S_maps / SO4_maps
final_comparison_both 

ggsave("Plots/final_comparison_env_maps.png", width = 15, height = 19, dpi = 600)



######################
# Category line charts
######################

epme_data_pre$stg <- factor(epme_data_pre$stg, levels = c("48", "49", "50", "51", "52", "53", "54"))
epme_data_post$stg <- factor(epme_data_post$stg, levels = c("48", "49", "50", "51", "52", "53", "54"))

stage_order <- c("Wordian", "Capitanian", "Wuchiapingian", "Changhsingian", "Induan", "Olenekian", "Anisian")

# Phylum 
percentage_data <- epme_data %>%
  group_by(stage, phylum) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  group_by(stage) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ungroup()

percentage_data$stage <- factor(percentage_data$stage, levels = stage_order)

phylum_line_plot <- ggplot(percentage_data, aes(x = stage, 
                                                y = percent, 
                                                color = phylum, 
                                                group = phylum)) +
  geom_line(size = 1.8, alpha = 0.8) +
  geom_point(size = 3) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 4.6, y = 90, label = "PT Boundary", color = "red", hjust = 0, alpha = 0.6) +
  scale_color_viridis_d(option = "turbo") + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    
  ) +
  labs(
    title = "A)",
    x = "Geological Stage",
    y = "Relative Abundance (%)",
    color = "Phylum"
  )

# Motility
percentage_data <- epme_data %>%
  group_by(stage, motility) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  group_by(stage) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ungroup()

percentage_data$stage <- factor(percentage_data$stage, levels = stage_order)


motility_line_plot <- ggplot(percentage_data, aes(x = stage, 
                                                  y = percent, 
                                                  color = motility, 
                                                  group = motility)) +
  geom_line(size = 1.8, alpha = 0.8) +
  geom_point(size = 3) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 4.7, y = 90, label = "PT Boundary", color = "red", hjust = 0, alpha = 0.6) +
  scale_color_viridis_d(option = "turbo") + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    
  ) +
  labs(
    title = "B)",
    y = "Relative Abundance (%)",
    x = "Geological Stage",
    color = "Motility"
  )

# life habit
percentage_data <- epme_data %>%
  group_by(stage, life_habit) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  group_by(stage) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ungroup()

percentage_data$stage <- factor(percentage_data$stage, levels = stage_order)


life_habit_line_plot <- ggplot(percentage_data, aes(x = stage, 
                                                    y = percent, 
                                                    color = life_habit, 
                                                    group = life_habit)) +
  geom_line(size = 1.8, alpha = 0.8) +
  geom_point(size = 3) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 4.7, y = 90, label = "PT Boundary", color = "red", hjust = 0, alpha = 0.6) +
  scale_color_viridis_d(option = "turbo") + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    
  ) +
  labs(
    title = "E)",
    y = "Relative Abundance (%)",
    x = "Geological Stage",
    color = "Life habit"
  )

# Vision
percentage_data <- epme_data %>%
  filter(vision != "unknown") %>%
  group_by(stage, vision) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  group_by(stage) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ungroup()

percentage_data$stage <- factor(percentage_data$stage, levels = stage_order)


vision_line_plot <- ggplot(percentage_data, aes(x = stage, 
                                                y = percent, 
                                                color = vision, 
                                                group = vision)) +
  geom_line(size = 1.8, alpha = 0.8) +
  geom_point(size = 3) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 4.7, y = 90, label = "PT Boundary", color = "red", hjust = 0, alpha = 0.6) +
  scale_color_viridis_d(option = "turbo") + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    
  ) +
  labs(
    title = "C)",
    y = "Relative Abundance (%)",
    x = "Geological Stage",
    color = "Vision"
  )

# diet
percentage_data <- epme_data %>%
  group_by(stage, diet) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  group_by(stage) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  ungroup()

percentage_data$stage <- factor(percentage_data$stage, levels = stage_order)

diet_line_plot <- ggplot(percentage_data, aes(x = stage, 
                                              y = percent, 
                                              color = diet, 
                                              group = diet)) +
  geom_line(size = 1.8, alpha = 0.8) +
  geom_point(size = 3) +
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 4.7, y = 90, label = "PT Boundary", color = "red", hjust = 0, alpha = 0.6) +
  scale_color_viridis_d(option = "turbo") + 
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    
  ) +
  labs(
    title = "D)",
    y = "Relative Abundance (%)",
    x = "Geological Stage",
    color = "Diet"
  )

design <- "
  AA
  BC
  DE
"

# Assemble plots
line_plots_final <- phylum_line_plot + 
  motility_line_plot + vision_line_plot + 
  diet_line_plot + life_habit_line_plot + 
  plot_layout(design = design)

ggsave("Plots/line_plots_final.png", width = 16, height = 20, dpi = 600)

#########################################################################################################
# Plots - results
#########################################################################################################


##########################
# Box plot - pre vs post 
##########################

my_colors <- c(
  "temp" = "gold",
  "O2"   = "cadetblue",
  "pH"   = "red3",
  "H2S"  = "darkorange1",
  "POC"  = "green4",
  "SO4"  = "purple"
)
# Prep Data
pre_data <- ENM_results_pre %>% filter(test_AUC >= 0.8 & test_AUC < 1) %>%
  mutate(Epoch = "Pre-Extinction") %>% dplyr::select(genus, Epoch, ends_with(".permutation.importance"))

post_data <- ENM_results_post %>% filter(test_AUC >= 0.8 & test_AUC < 1) %>%
  mutate(Epoch = "Post-Extinction") %>% dplyr::select(genus, Epoch, ends_with(".permutation.importance"))

############ Seperate plots
##Pre

pre_data_long <- pre_data %>% pivot_longer(
  cols = starts_with("temp"):ends_with(".importance"), 
  names_to = "Variable",
  values_to = "Importance"
) %>%
  mutate(Variable = gsub(".permutation.importance", "", Variable))

ggplot(pre_data_long, aes(x = reorder(Variable, Importance, FUN = median), y = Importance, fill = Variable)) +
  geom_boxplot(
    alpha = 0.7, 
    outlier.shape = NA
  ) +
  scale_x_discrete(labels = c(
    "temp" = "Temperature",
    "POC"  = "Particulate Organic Carbon",
    "SO4"  = "Sulfate",
    "pH"   = "pH",
    "O2"   = "Diss. Oxygen",
    "H2S"  = "Hydrogen Sulfide"
  )) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 50, 100)) +
  coord_flip() +
  labs(
    title = "ENM results pre",
    subtitle = NULL,
    x = "Environmental Variable",
    y = "Permutation Importance (%)"
  ) +
  theme_bw()

############ Seperate plots
##Post

post_data_long <- post_data %>% pivot_longer(
  cols = starts_with("temp"):ends_with(".importance"), 
  names_to = "Variable",
  values_to = "Importance"
) %>%
  mutate(Variable = gsub(".permutation.importance", "", Variable))

ggplot(post_data_long, aes(x = reorder(Variable, Importance, FUN = median), y = Importance, fill = Variable)) +
  geom_boxplot(
    alpha = 0.7, 
    outlier.shape = NA
  ) +
  scale_x_discrete(labels = c(
    "temp" = "Temperature",
    "POC"  = "Particulate Organic Carbon",
    "SO4"  = "Sulfate",
    "pH"   = "pH",
    "O2"   = "Diss. Oxygen",
    "H2S"  = "Hydrogen Sulfide"
  )) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 50, 100)) +
  coord_flip() +
  labs(
    title = "ENM results post",
    subtitle = NULL,
    x = "Environmental Variable",
    y = "Permutation Importance (%)"
  ) +
  theme_bw()


# Combine
combined_importance <- bind_rows(pre_data, post_data) %>%
  pivot_longer(
    cols = starts_with("temp"):ends_with(".importance"), 
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "", Variable))

combined_importance$Epoch <- factor(combined_importance$Epoch, levels = c("Post-Extinction", "Pre-Extinction"))

ggplot(combined_importance, aes(x = reorder(Variable, Importance, FUN = median), y = Importance, fill = Variable)) +
  geom_boxplot_pattern(
    aes(pattern = Epoch),
    pattern_color = "black",     
    pattern_fill = "black",
    pattern_angle = 45,          
    pattern_density = 0.2,       
    pattern_spacing = 0.03,      
    alpha = 0.7, 
    outlier.shape = NA
  ) +
  scale_pattern_manual(values = c("Pre-Extinction" = "none", "Post-Extinction" = "stripe")) +
  scale_x_discrete(labels = c(
    "temp" = "Temperature",
    "DOM"  = "Diss. Organic Matter",
    "SO4"  = "Sulfate",
    "pH"   = "pH",
    "O2"   = "Diss. Oxygen",
    "H2S"  = "Hydrogen Sulfide"
  )) +
    scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  guides(
    fill = "none", 
    pattern = guide_legend(
      title = "Age",
      override.aes = list(
        fill = "white",      
        pattern = c("stripe", "none"),
        alpha = 1           
      )
    )
  ) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 50, 100)) +
  coord_flip() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Environmental Variable",
    y = "Permutation Importance (%)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "white"),
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(face = "bold")
        )

ggsave(filename = "Environmental drivers of distribution Pre vs Post Extinction.png",path =  "C:/Dissertation/Dissertation/Actual final plots", width = 10, height = 8)

#########################
# Without colinearity 
#########################

# No temp (*pre*)
ENM_results_pre_sans_temp_filtered_auc <- ENM_results_pre_sans_temp %>%
  filter(test_AUC >= 0.8 & test_AUC < 1)

# Combine them back together
importance_long_pre_sans_temp <- ENM_results_pre_sans_temp_filtered_auc %>%
  dplyr::select(genus, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(genus),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable)) 

# Skip if no data
if (nrow(importance_long_pre_sans_temp) == 0) {
  next
}

comparisons_o2 <- list(
  c("O2", "H2S"),
  c("O2", "pH"), 
  c("O2", "SO4"),
  c("O2", "DOM")
)

my_colors <- c(
  "temp" = "#F5D94B",
  "O2"   = "cadetblue",
  "pH"   = "#CC1A37",
  "H2S"  = "#FF7E00",
  "DOM"  = "green4",
  "SO4"  = "#5F006E"
)

pre_sans_temp_plot <- ggplot(importance_long_pre_sans_temp, aes(x = reorder(Variable, Importance, FUN = median),
                                          y = Importance, 
                                          fill = Variable)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_jitter(
    alpha = 0.3, 
    size = 1,
    position = position_jitterdodge(dodge.width = 0.8)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(limits = c("pH", "H2S", "SO4", "DOM", "O2", "temp")) +
  scale_x_discrete(labels = c("O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_y_continuous(
    trans = "pseudo_log", 
    breaks = c(0, 1, 10, 50, 100),
    labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Environmental Variable",
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(face = "bold"))


ggsave("Environmental drivers pre sans temp.png", path =  "C:/Dissertation/Dissertation")

niche_summary_sans_temp <- importance_long_pre_sans_temp %>%
  group_by(Variable) %>%
  summarise(
    Mean = mean(Importance, na.rm = TRUE),
    Median = median(Importance, na.rm = TRUE),
    SD = sd(Importance, na.rm = TRUE),
    Min = min(Importance, na.rm = TRUE),
    Max = max(Importance, na.rm = TRUE),
    Breadth_IQR = IQR(Importance, na.rm = TRUE), 
    Sample_Size = n()
  ) %>%
  ungroup()


# No o2 (*pre*)

comparisons_temp <- list(
  c("temp", "H2S"),
  c("temp", "pH"), 
  c("temp", "SO4"),
  c("temp", "DOM")
)

ENM_results_pre_sans_o2_filtered_auc <- ENM_results_pre_sans_o2 %>%
  filter(test_AUC >= 0.8 & test_AUC < 1)

# Combine them 
importance_long_pre_sans_o2 <- ENM_results_pre_sans_o2_filtered_auc %>%
  dplyr::select(genus, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(genus),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable)) 


pre_sans_o2_plot <- ggplot(importance_long_pre_sans_o2, aes(x = reorder(Variable, Importance, FUN = median),
                                        y = Importance, 
                                        fill = Variable)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_jitter(
    alpha = 0.3, 
    size = 1,
    position = position_jitterdodge(dodge.width = 0.8)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(limits = c("pH", "H2S", "SO4", "DOM", "O2", "temp")) +
  scale_x_discrete(labels = c("temp" = "Temperature",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_y_continuous(
    trans = "pseudo_log", 
    breaks = c(0, 1, 10, 50, 100),
    labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(face = "bold"))

ggsave("Environmental drivers pre sans o2.png", path =  "C:/Dissertation/Dissertation")

# Without o2 (*post*)

comparisons_temp2 <- list(
  c("temp", "H2S"),
  c("temp", "pH"), 
  c("temp", "SO4"),
  c("temp", "DOM")
)

ENM_results_post_sans_o2_filtered_auc <- ENM_results_post_sans_o2 %>%
  filter(test_AUC >= 0.8 & test_AUC < 1)

# Combine them
importance_long_post_sans_o2 <- ENM_results_post_sans_o2_filtered_auc %>%
  dplyr::select(genus, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(genus),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable)) 


post_sans_o2_plot <- ggplot(importance_long_post_sans_o2, aes(x = reorder(Variable, Importance, FUN = median),
                                        y = Importance, 
                                        fill = Variable)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_jitter(
    alpha = 0.3, 
    size = 1,
    position = position_jitterdodge(dodge.width = 0.8)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(limits = c("pH", "H2S", "SO4", "DOM", "O2", "temp")) +
  scale_x_discrete(labels = c("temp" = "Temperature",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_y_continuous(
    trans = "pseudo_log", 
    breaks = c(0, 1, 10, 50, 100),
    labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Importance (%)"
  ) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(face = "bold"))

ggsave("Drivers post extinction without o2.png", width = 14, height = 8)


# Without pH (*post*)

comparisons_temp <- list(
  c("temp", "H2S"),
  c("temp", "O2"), 
  c("temp", "SO4"),
  c("temp", "DOM")
)

ENM_results_post_sans_ph_filtered_auc <- ENM_results_post_sans_ph %>%
  filter(test_AUC >= 0.8 & test_AUC < 1)

# Combine them
importance_long_post_sans_ph <- ENM_results_post_sans_ph_filtered_auc %>%
  dplyr::select(genus, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(genus),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable)) 


post_sans_ph_plot <- ggplot(importance_long_post_sans_ph, aes(x = reorder(Variable, Importance, FUN = median),
                                         y = Importance, 
                                         fill = Variable)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_jitter(
    alpha = 0.3, 
    size = 1,
    position = position_jitterdodge(dodge.width = 0.8)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_x_discrete(limits = c("pH", "H2S", "SO4", "DOM", "O2", "temp")) +
  scale_x_discrete(labels = c("temp" = "Temperature",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_y_continuous(
    trans = "pseudo_log", 
    breaks = c(0, 1, 10, 50, 100),
    labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Environmental Variable",
    y = "Importance (%)"
  ) +
  theme_bw() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(face = "bold"))

ggsave("Drivers post extinction without pH.png", width = 14, height = 8)

top_row <- (pre_sans_temp_plot + labs(tag = "a")) + (pre_sans_o2_plot + labs(tag = "b")) + plot_layout(ncol = 2)
bot_row <- (post_sans_ph_plot + labs(tag = "c")) + (post_sans_o2_plot + labs(tag = "c")) + plot_layout(ncol = 2)

ENV_vars_boxplots_all <- (
  wrap_elements(full = top_row + plot_annotation(title = "Pre-extinction", theme = theme(plot.title = element_text(hjust = 0.52, face = "bold")))) /
    wrap_elements(full = bot_row + plot_annotation(title = "Post-extinction", theme = theme(plot.title = element_text(hjust = 0.52, face = "bold")))) 
  )


ENV_vars_boxplots_all
ggsave("ENV_vars_boxplots_all.png", width = 14, height = 8)


#########################
# Correlation grids
#########################

cor_pre <- epme_data_pre_colintest %>% select(temp, O2, pH, H2S, DOM, SO4, sal, ca, ALK, DIC, Fe, Mg, L, FeL)
cor1 = cor(cor_pre, method = "spearman")
corrplot.mixed(cor1, lower.col = "black", number.cex = .7)

cor_post <- epme_data_post_colintest_post %>% select(temp, O2, pH, H2S, DOM, SO4, sal, ca, ALK, DIC, Fe, Mg, L, FeL)
cor2 = cor(cor_post, method = "spearman")
corrplot.mixed(cor2, lower.col = "black", number.cex = .7)

#######################################
# Plots separated by categories
#######################################

phylum_high_auc <- ENM_results_pre %>%
  filter(test_AUC >= 0.8) %>%
  filter(phylum %in% c("Chordata", "Brachiopoda", "Porifera", "Mollusca", 
                       "Annelida", "Foraminifera", "Cnidaria", 
                       "Echinodermata", "Bryozoa", "Arthropoda"))



plot_data <- phylum_high_auc %>%
  dplyr::select(phylum, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(phylum),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

phylum_comparison_plot <- ggplot(plot_data, aes(x = Importance, y = Variable, fill = phylum)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Pre-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = "Environmental Variable",
    fill = "Phylum"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )


phylum_high_auc_post <- ENM_results_post %>%
  filter(test_AUC >= 0.8) %>%
  filter(phylum %in% c("Chordata", "Brachiopoda", "Porifera", "Mollusca", 
                       "Annelida", "Foraminifera", "Cnidaria", 
                       "Echinodermata", "Bryozoa", "Arthropoda"))



plot_data_post <- phylum_high_auc_post %>%
  dplyr::select(phylum, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(phylum),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

phylum_comparison_plot_post <- ggplot(plot_data_post, aes(x = Importance, y = Variable, fill = phylum)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Post-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = NULL,
    fill = "Phylum"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

phylum_plots <- phylum_comparison_plot + phylum_comparison_plot_post


# Motility
motility_high_auc <- ENM_results_pre %>%
  filter(test_AUC >= 0.8) %>%
  filter(motility %in% c("stationary, attached", "stationary", 
                         "facultatively mobile", "slow-moving", 
                         "stationary, attached, epibiont", "fast-moving", 
                         "actively mobile"))

plot_data_pre <- motility_high_auc %>%
  dplyr::select(motility, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(motility),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

motility_comparison_plot <- ggplot(plot_data_pre, aes(x = Importance, y = Variable, fill = motility)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Pre-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = "Environmental Variable",
    fill = "motility"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

motility_high_auc_post <- ENM_results_post %>%
  filter(test_AUC >= 0.8) %>%
  filter(motility %in% c("stationary, attached", "stationary", 
                         "facultatively mobile", "slow-moving", 
                         "stationary, attached, epibiont", "fast-moving", 
                         "actively mobile"))

plot_data_post <- motility_high_auc_post %>%
  dplyr::select(motility, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(motility),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

motility_comparison_plot_post <- ggplot(plot_data_post, aes(x = Importance, y = Variable, fill = motility)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Post-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = NULL,
    fill = "motility"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

motility_plots <- motility_comparison_plot + motility_comparison_plot_post

####################
# some stats tests
####################

actively_mobile <- plot_data_pre %>% filter(motility == "actively mobile")
facultatively_mobile <- plot_data_pre %>% filter(motility == "facultatively mobile")
stationary <-  plot_data_pre %>% filter(motility == "stationary")


temp_data_pre <- plot_data_pre %>% filter(Variable == "temp")

wilcox.test(Importance ~ motility, 
            data = temp_data_pre %>% filter(motility %in% c("actively mobile", "facultatively mobile")))

wilcox.test(stationary$Variable)

act_dunn<-dunnTest(Importance ~ Variable, data = actively_mobile)
facul_dunn<-dunnTest(Importance ~ Variable, data = facultatively_mobile)
stat_dunn<-dunnTest(Importance ~ Variable, data = stationary)


niche_test_results <- combined_niche %>%
  group_by(Variable) %>%
   do(tidy(wilcox.test(Value ~ Period, data = .)))


active_post <- plot_data_post %>% filter(motility == "actively mobile")
facul_post <- plot_data_post %>% filter(motility == "facultatively mobile")
stationary_post <- plot_data_post %>% filter(motility == "stationary")

act_dunn_post<-dunnTest(Importance ~ Variable, data = active_post)
facul_dunn_post<-dunnTest(Importance ~ Variable, data = facul_post)
stat_dunn_post<-dunnTest(Importance ~ Variable, data = stationary_post)

# Vision
vision_high_auc <- ENM_results_pre %>%
  filter(test_AUC >= 0.8) %>%
  filter(vision %in% c("blind", "limited", "well-developed"))

plot_data <- vision_high_auc %>%
  dplyr::select(vision, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(vision),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

vision_comparison_plot <- ggplot(plot_data, aes(x = Importance, y = Variable, fill = vision)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Pre-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = "Environmental Variable",
    fill = "vision"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

vision_high_auc_post <- ENM_results_post %>%
  filter(test_AUC >= 0.8) %>%
  filter(vision %in% c("blind", "limited", "well-developed"))

plot_data_post <- vision_high_auc_post %>%
  dplyr::select(vision, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(vision),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

vision_comparison_plot_post <- ggplot(plot_data_post, aes(x = Importance, y = Variable, fill = vision)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Post-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = NULL,
    fill = "vision"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

vision_plots <- vision_comparison_plot + vision_comparison_plot_post

blind <- plot_data %>% filter(vision == "blind")
limited <- plot_data %>% filter(vision == "limited")
well_developed <- plot_data %>% filter(vision == "well-developed")
unknown <- plot_data %>% filter(vision == "unknown")

#############
# Stats tests
#############

blind_dunn<-dunnTest(Importance ~ Variable, data = blind)
limited_dunn<-dunnTest(Importance ~ Variable, data = limited)
well_developed_dunn<-dunnTest(Importance ~ Variable, data = well_developed)
unknown_dunn<-dunnTest(Importance ~ Variable, data = unknown)

blind_post <- plot_data_post %>% filter(vision == "blind")
limited_post <- plot_data_post %>% filter(vision == "limited")
well_developed_post <- plot_data_post %>% filter(vision == "well-developed")
unknown_post <- plot_data_post %>% filter(vision == "unknown")

blind_dunn_post<-dunnTest(Importance ~ Variable, data = blind_post)
limited_dunn_post<-dunnTest(Importance ~ Variable, data = limited_post)
well_developed_dunn_post<-dunnTest(Importance ~ Variable, data = well_developed_post)
unknown_dunn_post<-dunnTest(Importance ~ Variable, data = unknown_post)

blind_enm <- vision_high_auc %>% filter(vision == "blind")
blind_enm_post <- vision_high_auc_post %>% filter(vision == "blind")
limited_enm <- vision_high_auc %>% filter(vision == "limited")
limited_enm_post <- vision_high_auc_post %>% filter(vision == "limited")
well_developed_enm <- vision_high_auc %>% filter(vision == "well-developed")
well_developed_enm_post <- vision_high_auc_post %>% filter(vision == "well-developed")
unknown_enm <- vision_high_auc %>% filter(vision == "unknown")
unknown_enm_post <- vision_high_auc_post %>% filter(vision == "unknown")


wilcox.test(blind_enm$temp.permutation.importance, blind_enm_post$temp.permutation.importance)
wilcox.test(blind_enm$O2.permutation.importance, blind_enm_post$O2.permutation.importance)
wilcox.test(blind_enm$DOM.permutation.importance, blind_enm_post$DOM.permutation.importance)
wilcox.test(blind_enm$pH.permutation.importance, blind_enm_post$pH.permutation.importance)
wilcox.test(blind_enm$SO4.permutation.importance, blind_enm_post$SO4.permutation.importance)
wilcox.test(blind_enm$H2S.permutation.importance, blind_enm_post$H2S.permutation.importance)

wilcox.exact(limited_enm$temp.permutation.importance, limited_enm_post$temp.permutation.importance)
wilcox.exact(limited_enm$O2.permutation.importance, limited_enm_post$O2.permutation.importance)
wilcox.exact(limited_enm$DOM.permutation.importance, limited_enm_post$DOM.permutation.importance)
wilcox.exact(limited_enm$pH.permutation.importance, limited_enm_post$pH.permutation.importance)
wilcox.exact(limited_enm$SO4.permutation.importance, limited_enm_post$SO4.permutation.importance)
wilcox.exact(limited_enm$H2S.permutation.importance, limited_enm_post$H2S.permutation.importance)

wilcox.test(well_developed_enm$temp.permutation.importance, well_developed_enm_post$temp.permutation.importance)
wilcox.test(well_developed_enm$O2.permutation.importance, well_developed_enm_post$O2.permutation.importance)
wilcox.test(well_developed_enm$DOM.permutation.importance, well_developed_enm_post$DOM.permutation.importance)
wilcox.test(well_developed_enm$SO4.permutation.importance, well_developed_enm_post$SO4.permutation.importance)
wilcox.test(well_developed_enm$H2S.permutation.importance, well_developed_enm_post$H2S.permutation.importance)
wilcox.test(well_developed_enm$pH.permutation.importance, well_developed_enm_post$pH.permutation.importance)

wilcox.exact(unknown_enm$temp.permutation.importance, unknown_enm_post$temp.permutation.importance)
wilcox.exact(unknown_enm$O2.permutation.importance, unknown_enm_post$O2.permutation.importance)
wilcox.exact(unknown_enm$DOM.permutation.importance, unknown_enm_post$DOM.permutation.importance)
wilcox.exact(unknown_enm$pH.permutation.importance, unknown_enm_post$pH.permutation.importance)
wilcox.exact(unknown_enm$SO4.permutation.importance, unknown_enm_post$SO4.permutation.importance)
wilcox.exact(unknown_enm$H2S.permutation.importance, unknown_enm_post$H2S.permutation.importance)


# Diet
diet_high_auc <- ENM_results_pre %>%
  filter(test_AUC >= 0.8) %>%
  filter(diet %in% c("suspension feeder", "other", "carnivore", "omnivore", "grazer" ))

plot_data_pre <- diet_high_auc %>%
  dplyr::select(diet, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(diet),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

diet_comparison_plot <- ggplot(plot_data_pre, aes(x = Importance, y = Variable, fill = diet)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Pre-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = "Environmental Variable",
    fill = "diet"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

diet_high_auc_post <- ENM_results_post %>%
  filter(test_AUC >= 0.8) %>%
  filter(diet %in% c("suspension feeder", "other", "carnivore", "omnivore", "grazer"))

plot_data_post <- diet_high_auc_post %>%
  dplyr::select(diet, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(diet),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

diet_comparison_plot_post <- ggplot(plot_data_post, aes(x = Importance, y = Variable, fill = diet)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Post-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = NULL,
    fill = "diet"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

diet_plots <- diet_comparison_plot + diet_comparison_plot_post

carnivore <- plot_data_pre %>% filter(diet == "carnivore")
omnivore <- plot_data_pre %>% filter(diet == "omnivore")
suspension_feeder <- plot_data_pre %>% filter(diet == "suspension feeder")
other <- plot_data_pre %>% filter(diet == "other")

#############
# Stats tests
#############

carnivore_dunn<-dunnTest(Importance ~ Variable, data = carnivore)
omnivore_dunn<-dunnTest(Importance ~ Variable, data = omnivore)
suspension_feeder_dunn<-dunnTest(Importance ~ Variable, data = suspension_feeder)
other_dunn<-dunnTest(Importance ~ Variable, data = other)

carnivore_post <- plot_data_post %>% filter(diet == "carnivore")
omnivore_post <- plot_data_post %>% filter(diet == "omnivore")
suspension_feeder_post <- plot_data_post %>% filter(diet == "suspension feeder")
other_post <- plot_data_post %>% filter(diet == "other")

carnivore_dunn_post<-dunnTest(Importance ~ Variable, data = carnivore_post)
omnivore_dunn_post<-dunnTest(Importance ~ Variable, data = omnivore_post)
suspension_feeder_dunn_post<-dunnTest(Importance ~ Variable, data = suspension_feeder_post)
other_dunn_post<-dunnTest(Importance ~ Variable, data = other_post)

wilcox.test(blind$Importance, blind_post$Importance)
wilcox.test(limited$Importance, limited_post$Importance)
wilcox.test(well_developed$Importance, well_developed_post$Importance)


carnivore_enm <- diet_high_auc %>% filter(diet == "carnivore")
carnivore_enm_post <- diet_high_auc_post %>% filter(diet == "carnivore")
omnivore_enm <- diet_high_auc %>% filter(diet == "omnivore")
omnivore_enm_post <- diet_high_auc_post %>% filter(diet == "omnivore")
suspension_feeder_enm <- diet_high_auc %>% filter(diet == "suspension feeder")
suspension_feeder_enm_post <- diet_high_auc_post %>% filter(diet == "suspension feeder")
other_enm <- diet_high_auc %>% filter(diet == "other")
other_enm_post <- diet_high_auc_post %>% filter(diet == "other")


wilcox.test(carnivore_enm$temp.permutation.importance, carnivore_enm_post$temp.permutation.importance)
wilcox.test(carnivore_enm$O2.permutation.importance, carnivore_enm_post$O2.permutation.importance)
wilcox.test(carnivore_enm$DOM.permutation.importance, carnivore_enm_post$DOM.permutation.importance)
wilcox.test(carnivore_enm$pH.permutation.importance, carnivore_enm_post$pH.permutation.importance)
wilcox.test(carnivore_enm$SO4.permutation.importance, carnivore_enm_post$SO4.permutation.importance)
wilcox.test(carnivore_enm$H2S.permutation.importance, carnivore_enm_post$H2S.permutation.importance)

wilcox.exact(omnivore_enm$temp.permutation.importance, omnivore_enm_post$temp.permutation.importance)
wilcox.exact(omnivore_enm$O2.permutation.importance, omnivore_enm_post$O2.permutation.importance)
wilcox.exact(omnivore_enm$DOM.permutation.importance, omnivore_enm_post$DOM.permutation.importance)
wilcox.exact(omnivore_enm$pH.permutation.importance, omnivore_enm_post$pH.permutation.importance)
wilcox.exact(omnivore_enm$SO4.permutation.importance, omnivore_enm_post$SO4.permutation.importance)
wilcox.exact(omnivore_enm$H2S.permutation.importance, omnivore_enm_post$H2S.permutation.importance)

wilcox.test(suspension_feeder_enm$temp.permutation.importance, suspension_feeder_enm_post$temp.permutation.importance)
wilcox.test(suspension_feeder_enm$O2.permutation.importance, suspension_feeder_enm_post$O2.permutation.importance)
wilcox.test(suspension_feeder_enm$DOM.permutation.importance, suspension_feeder_enm_post$DOM.permutation.importance)
wilcox.test(suspension_feeder_enm$SO4.permutation.importance, suspension_feeder_enm_post$SO4.permutation.importance)
wilcox.test(suspension_feeder_enm$H2S.permutation.importance, suspension_feeder_enm_post$H2S.permutation.importance)
wilcox.test(suspension_feeder_enm$pH.permutation.importance, suspension_feeder_enm_post$pH.permutation.importance)

wilcox.exact(other_enm$temp.permutation.importance, other_enm_post$temp.permutation.importance)
wilcox.exact(other_enm$O2.permutation.importance, other_enm_post$O2.permutation.importance)
wilcox.exact(other_enm$DOM.permutation.importance, other_enm_post$DOM.permutation.importance)
wilcox.exact(other_enm$pH.permutation.importance, other_enm_post$pH.permutation.importance)
wilcox.exact(other_enm$SO4.permutation.importance, other_enm_post$SO4.permutation.importance)
wilcox.exact(other_enm$H2S.permutation.importance, other_enm_post$H2S.permutation.importance)



# Life habit
life_habit_high_auc <- ENM_results_pre %>%
  filter(test_AUC >= 0.8) %>%
  filter(life_habit %in% c("epifaunal", "unknown", "infaunal", "nektonic", "nektobenthic"))

plot_data_pre <- life_habit_high_auc %>%
  dplyr::select(life_habit, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(life_habit),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

life_habit_comparison_plot <- ggplot(plot_data_pre, aes(x = Importance, y = Variable, fill = life_habit)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Pre-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = "Environmental Variable",
    fill = "life_habit"
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

life_habit_high_auc_post <- ENM_results_post %>%
  filter(test_AUC >= 0.8) %>%
  filter(life_habit %in% c("epifaunal", "unknown", "infaunal", "nektonic", "nektobenthic"))

plot_data_post <- life_habit_high_auc_post %>%
  dplyr::select(life_habit, ends_with(".permutation.importance")) %>%
  pivot_longer(
    cols = -c(life_habit),
    names_to = "Variable",
    values_to = "Importance"
  ) %>%
  mutate(Variable = gsub(".permutation.importance", "",Variable))

life_habit_comparison_plot_post <- ggplot(plot_data_post, aes(x = Importance, y = Variable, fill = life_habit)) +
  geom_vline(xintercept = c(0, 1, 10, 50, 100), color = "grey90") +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.8) +
  scale_y_discrete(labels = c("temp" = "Temperature",
                              "O2" = "Diss. oxygen",
                              "DOM" = "Diss. organic matter",
                              "H2S" = "Hydrogen sulphide",
                              "SO4" = "Sulphate")) +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 10, 50, 100)) +
  theme_bw() +
  labs(
    title = "Post-extinction",
    subtitle = NULL,
    x = "Permutation Importance (%)",
    y = NULL,
    fill = "life habit"
  ) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

life_habit_plots <- life_habit_comparison_plot + life_habit_comparison_plot_post

#############
# Stats tests
#############

epifaunal <- plot_data_pre %>% filter(life_habit == "epifaunal")
infaunal <- plot_data_pre %>% filter(life_habit == "infaunal")
nektobenthic <- plot_data_pre %>% filter(life_habit == "nektobenthic")
nektonic <- plot_data_pre %>% filter(life_habit == "nektonic")
unknown <- plot_data_pre %>% filter(life_habit == "unknown")

epifaunal_dunn<-dunnTest(Importance ~ Variable, data = epifaunal)
infaunal_dunn<-dunnTest(Importance ~ Variable, data = infaunal)
nektobenthic_dunn<-dunnTest(Importance ~ Variable, data = nektobenthic)
nektonic_dunn<-dunnTest(Importance ~ Variable, data = nektonic)
unknown_dunn<-dunnTest(Importance ~ Variable, data = unknown)

epifaunal_post <- plot_data_post %>% filter(life_habit == "epifaunal")
infaunal_post <- plot_data_post %>% filter(life_habit == "infaunal")
nektobenthic_post <- plot_data_post %>% filter(life_habit == "nektobenthic")
nektonic_post <- plot_data_post %>% filter(life_habit == "nektonic")
unknown_post <- plot_data_post %>% filter(life_habit == "unknown")

epifaunal_dunn_post<-dunnTest(Importance ~ Variable, data = epifaunal_post)
infaunal_dunn_post<-dunnTest(Importance ~ Variable, data = infaunal_post)
nektobenthic_dunn_post<-dunnTest(Importance ~ Variable, data = nektobenthic_post)
nektonic_dunn_post<-dunnTest(Importance ~ Variable, data = nektonic_post)
unknown_dunn_post<-dunnTest(Importance ~ Variable, data = unknown_post)


epifaunal_enm <- life_habit_high_auc %>% filter(life_habit == "epifaunal")
epifaunal_enm_post <- life_habit_high_auc_post %>% filter(life_habit == "epifaunal")
infaunal_enm <- life_habit_high_auc %>% filter(life_habit == "infaunal")
infaunal_enm_post <- life_habit_high_auc_post %>% filter(life_habit == "infaunal")
nektobenthic_enm <- life_habit_high_auc %>% filter(life_habit == "nektobenthic")
nektobenthic_enm_post <- life_habit_high_auc_post %>% filter(life_habit == "nektobenthic")
nektonic_enm <- life_habit_high_auc %>% filter(life_habit == "nektonic")
nektonic_enm_post <- life_habit_high_auc_post %>% filter(life_habit == "nektonic")
unknown_enm <- life_habit_high_auc %>% filter(life_habit == "unknown")
unknown_enm_post <- life_habit_high_auc_post %>% filter(life_habit == "unknown")



wilcox.test(epifaunal_enm$temp.permutation.importance, epifaunal_enm_post$temp.permutation.importance)
wilcox.test(epifaunal_enm$O2.permutation.importance, epifaunal_enm_post$O2.permutation.importance)
wilcox.test(epifaunal_enm$DOM.permutation.importance, epifaunal_enm_post$DOM.permutation.importance)
wilcox.test(epifaunal_enm$pH.permutation.importance, epifaunal_enm_post$pH.permutation.importance)
wilcox.test(epifaunal_enm$SO4.permutation.importance, epifaunal_enm_post$SO4.permutation.importance)
wilcox.test(epifaunal_enm$H2S.permutation.importance, epifaunal_enm_post$H2S.permutation.importance)

wilcox.exact(infaunal_enm$temp.permutation.importance, infaunal_enm_post$temp.permutation.importance)
wilcox.exact(infaunal_enm$O2.permutation.importance, infaunal_enm_post$O2.permutation.importance)
wilcox.exact(infaunal_enm$DOM.permutation.importance, infaunal_enm_post$DOM.permutation.importance)
wilcox.exact(infaunal_enm$pH.permutation.importance, infaunal_enm_post$pH.permutation.importance)
wilcox.exact(infaunal_enm$SO4.permutation.importance, infaunal_enm_post$SO4.permutation.importance)
wilcox.exact(infaunal_enm$H2S.permutation.importance, infaunal_enm_post$H2S.permutation.importance)

wilcox.exact(nektobenthic_enm$temp.permutation.importance, nektobenthic_enm_post$temp.permutation.importance)
wilcox.exact(nektobenthic_enm$O2.permutation.importance, nektobenthic_enm_post$O2.permutation.importance)
wilcox.exact(nektobenthic_enm$DOM.permutation.importance, nektobenthic_enm_post$DOM.permutation.importance)
wilcox.exact(nektobenthic_enm$SO4.permutation.importance, nektobenthic_enm_post$SO4.permutation.importance)
wilcox.exact(nektobenthic_enm$H2S.permutation.importance, nektobenthic_enm_post$H2S.permutation.importance)
wilcox.exact(nektobenthic_enm$pH.permutation.importance, nektobenthic_enm_post$pH.permutation.importance)

wilcox.test(nektonic_enm$temp.permutation.importance, nektonic_enm_post$temp.permutation.importance)
wilcox.test(nektonic_enm$O2.permutation.importance, nektonic_enm_post$O2.permutation.importance)
wilcox.test(nektonic_enm$DOM.permutation.importance, nektonic_enm_post$DOM.permutation.importance)
wilcox.test(nektonic_enm$pH.permutation.importance, nektonic_enm_post$pH.permutation.importance)
wilcox.test(nektonic_enm$SO4.permutation.importance, nektonic_enm_post$SO4.permutation.importance)
wilcox.test(nektonic_enm$H2S.permutation.importance, nektonic_enm_post$H2S.permutation.importance)

wilcox.exact(unknown_enm$temp.permutation.importance, unknown_enm_post$temp.permutation.importance)
wilcox.exact(unknown_enm$O2.permutation.importance, unknown_enm_post$O2.permutation.importance)
wilcox.exact(unknown_enm$DOM.permutation.importance, unknown_enm_post$DOM.permutation.importance)
wilcox.exact(unknown_enm$pH.permutation.importance, unknown_enm_post$pH.permutation.importance)
wilcox.exact(unknown_enm$SO4.permutation.importance, unknown_enm_post$SO4.permutation.importance)
wilcox.exact(unknown_enm$H2S.permutation.importance, unknown_enm_post$H2S.permutation.importance)

##############################################
# paleolat graph
###############################################
# Ensure you have your full list of numeric stages for consistent axes
all_stages <- as.numeric(as.character(unique(global_survivor_epme_data$stg)))
phyla_list <- unique(global_survivor_epme_data$phylum)

for (p in phyla_list) {
  
  current_data <- global_survivor_epme_data %>% 
    filter(phylum == p) %>%
    mutate(stg_num = as.numeric(as.character(stg)))
  
  if(nrow(current_data) == 0) next
  
  # Lat plot
  lat_migration_plot <- ggplot(current_data, aes(x = stg_num, y = paleolat)) +
    geom_vline(xintercept = 51.5, linetype = "dashed", color = "red", linewidth = 1, alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_violin(aes(group = stg_num), alpha = 0.5, fill = "grey80", color = NA, scale = "width", bw = 2) +
    geom_jitter(alpha = 0.5, size = 1.5, color = "steelblue", width = 0.2) +
    annotate("text", x = 52.2, y = 87, label = "PT Boundary", 
             color = "red", alpha = 0.4, fontface = "bold", vjust = -0.5) +
    theme_minimal() +
    labs(
      title = paste("Latitudinal Distribution of", p, "Survivors"),
      subtitle = "Tracking migration away from the equator during the EPME ",
      x = "Geological Stage",
      y = "Paleolatitude (Degrees N/S)"
    ) +
    
    scale_x_continuous(
      breaks = all_stages, 
      limits = c(min(all_stages) - 0.5, max(all_stages) + 0.5)
    ) +
    
    scale_y_continuous(breaks = seq(-90, 90, 30), limits = c(-90, 90))
  
  # Lng plot
  lng_migration_plot <- ggplot(current_data, aes(x = stg_num, y = paleolng)) +
    geom_vline(xintercept = 51.5, linetype = "dashed", color = "red", linewidth = 1, alpha = 0.4) +
    geom_violin(aes(group = stg_num), alpha = 0.7, fill = "grey80", color = NA, scale = "width", bw = 2) +
    geom_jitter(alpha = 0.5, size = 1.5, color = "steelblue", width = 0.2) +
    coord_flip() +
    
    annotate("text", x = 51.5, y = 160, label = "PT Boundary", 
             color = "red", alpha = 0.4, fontface = "bold", vjust = -0.5) +
    
    theme_minimal() +
    labs(
      title = paste("Longitudinal Distribution of", p, "Survivors"),
      subtitle = "Tracking migration during the PTB crisis",
      x = "Geological Stage",
      y = "Paleolongitude (Degrees)"
    ) +
    
    scale_x_continuous(
      breaks = all_stages, 
      limits = c(min(all_stages) - 0.5, max(all_stages) + 0.5)
    ) +
    scale_y_continuous(breaks = seq(-180, 180, 60), limits = c(-180, 180))
  
  
  
  combined_plot <- lat_migration_plot + lng_migration_plot + 
    plot_annotation(
      title = paste("Spatial Distribution of", p, "Survivors"),
      subtitle = "Tracking migration across the Permian Triassic boundary",
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  file_name <- paste0("Combined_Spatial_Shift_", p, ".png")
  ggsave(file_name, combined_plot, width = 14, height = 7)
  
  message(paste("Successfully created plot for:", p))
}
################################
# Statistical tests
###############################
combined_importance_pre <- combined_importance %>% filter(Epoch == "Pre-Extinction")
combined_importance_post <- combined_importance %>% filter(Epoch == "Post-Extinction")

  dunn_output_pre <- dunnTest(Importance ~ Variable, data = combined_importance_pre)
  dunn_output_post <- dunnTest(Importance ~ Variable, data = combined_importance_post)
  
  aggregate(rank(combined_importance_pre$Importance) ~ Variable, 
            data = combined_importance_pre, 
            FUN = mean)
  
  wilcox.test(pre_data$temp.permutation.importance, post_data$temp.permutation.importance)
  wilcox.test(pre_data$DOM.permutation.importance, post_data$DOM.permutation.importance)
  wilcox.test(pre_data$SO4.permutation.importance, post_data$SO4.permutation.importance)
  wilcox.test(pre_data$pH.permutation.importance, post_data$pH.permutation.importance)
  wilcox.test(pre_data$O2.permutation.importance, post_data$O2.permutation.importance)
  wilcox.test(pre_data$H2S.permutation.importance, post_data$H2S.permutation.importance)
  
  
  ############
  dunn_output_pre_sans_temp <- dunnTest(Importance ~ Variable, data = importance_long_pre_sans_temp )
  dunn_output_pre_sans_o2 <- dunnTest(Importance ~ Variable, data = importance_long_pre_sans_o2 )
  dunn_output_post_sans_o2 <- dunnTest(Importance ~ Variable, data = importance_long_post_sans_o2 )
  dunn_output_post_sans_ph <- dunnTest(Importance ~ Variable, data = importance_long_post_sans_ph )
  
  ##################
  # Run the Wilcoxon Rank-Sum Test
  wilcox.test(ENM_results_pre_sans_o2_filtered_auc$temp.permutation.importance, ENM_results_post_sans_o2_filtered_auc$temp.permutation.importance)
  wilcox.test(ENM_results_pre_sans_o2_filtered_auc$DOM.permutation.importance, ENM_results_post_sans_o2_filtered_auc$DOM.permutation.importance)
  wilcox.test(ENM_results_pre_sans_o2_filtered_auc$SO4.permutation.importance, ENM_results_post_sans_o2_filtered_auc$SO4.permutation.importance)
  wilcox.test(ENM_results_pre_sans_o2_filtered_auc$H2S.permutation.importance, ENM_results_post_sans_o2_filtered_auc$H2S.permutation.importance)
  wilcox.test(ENM_results_pre_sans_o2_filtered_auc$pH.permutation.importance, ENM_results_post_sans_o2_filtered_auc$pH.permutation.importance)
  
  
  wilcox.test(ENM_results_pre_sans_temp_filtered_auc$O2.permutation.importance, ENM_results_post_sans_ph_filtered_auc$O2.permutation.importance)
  wilcox.test(ENM_results_pre_sans_temp_filtered_auc$DOM.permutation.importance, ENM_results_post_sans_ph_filtered_auc$DOM.permutation.importance)
  wilcox.test(ENM_results_pre_sans_temp_filtered_auc$SO4.permutation.importance, ENM_results_post_sans_ph_filtered_auc$SO4.permutation.importance)
  wilcox.test(ENM_results_pre_sans_temp_filtered_auc$H2S.permutation.importance, ENM_results_post_sans_ph_filtered_auc$H2S.permutation.importance)
  
  
  ##############################
  # Avg var data
  ##############################

  
  niche_data <- epme_data_pre %>%
    pivot_longer(cols = c(temp, O2, pH, DOM, SO4, H2S), 
                 names_to = "Variable", 
                 values_to = "Value")
  
  niche_data_post <- epme_data_post %>%
    pivot_longer(cols = c(temp, O2, pH, DOM, SO4, H2S), 
                 names_to = "Variable", 
                 values_to = "Value")

  
  phyla_list <- unique(niche_data$phylum)
  
  plot_list <- list()
  
  for (p in phyla_list) {
    
    temp_data <- niche_data %>% filter(phylum == p)
    
    p_plot <- ggplot(temp_data, aes(x = Value, fill = Variable)) +
      geom_histogram(aes(y = after_stat(density)), 
                     bins = 30,           
                     alpha = 0.8,         
                     position = "identity", 
                     color = "white") +
      geom_density(alpha = 0.5, colour = "black") +
      facet_wrap(~Variable, scales = "free", ncol = 1, strip.position = "left", axes = "all") +
      scale_fill_manual(values = my_colors) +
      labs(title = p, 
           x = NULL, 
           y = NULL) +
      theme_bw() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(size = 8)
      )
    
    if (p == phyla_list[1]) {
      p_plot <- p_plot + theme(
      strip.text.y = element_text(angle = 0, face = "bold", size = 10)
      )
    } else {
      p_plot <- p_plot + theme(
        strip.text.y = element_blank() 
      )
    }
    
    plot_list[[p]] <- p_plot
  }
  final_stitch <- plot_list[["Brachiopoda"]]|plot_list[["Mollusca"]]|plot_list[["Chordata"]]|plot_list[["Arthropoda"]]|plot_list[["Porifera"]]|plot_list[["Bryozoa"]]|plot_list[["Cnidaria"]]|plot_list[["Echinodermata"]]|plot_list[["Foraminifera"]]|plot_list[["Annelida"]]|plot_list[["Hyolitha"]]
  final_stitch
  
  #######################################################  
  # Some plots and stats
  ######################################################
  
  qq_vars <- c("temp", "H2S", "O2", "SO4", "pH", "DOM")
  
  qq_plot_list <- list()
  
  for (v in qq_vars) {
    
    p <- ggplot(epme_data_pre, aes(sample = .data[[v]])) +
      stat_qq(color = "blue", alpha = 0.5) + 
      stat_qq_line(color = "red") +
      labs(title = paste("Q-Q Plot:", v),
           x = "Theoretical",
           y = "Sample") +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"))
    
    qq_plot_list[[v]] <- p
  }
  
  combined_qq_plots <- wrap_plots(qq_plot_list, ncol = 3) + 
    plot_annotation(title = "Normality Assessment (Pre-extinction)",
                    subtitle = "Linear points = Normal; Curved points = Non-normal",
                    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  

  pres <- ggplot(niche_data, aes(x = Value, fill = Variable)) +
    geom_histogram(aes(y = after_stat(density)), 
                   bins = 30,           
                   alpha = 0.8,         
                   position = "identity", 
                   color = "white") +
    geom_density(alpha = 0.5, colour = "black") +
    facet_wrap(~Variable, scales = "free", nrow = 1) +
    labs(title = "Pre-extinction",
         x = NULL,
         y = "Occurrence Density") +
    theme_bw() +
    theme (legend.position = "none",
             plot.title = element_text(hjust = 0.5, face = "bold"),
             panel.grid.minor = element_blank()
             ) 
  
  posts <- ggplot(niche_data_post, aes(x = Value, fill = Variable)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Variable, scales = "free", nrow = 1) +
    labs(title = "Post-extinction",
         x = "Environmental Gradient (Scaled 0-100)",
         y = "Occurrence Density") +
    theme_bw() +
    theme (legend.position = "none",
           plot.title = element_text(hjust = 0.5, face = "bold"),
           panel.grid.minor = element_blank()
    )
  
  densitys <- pres / posts
  densitys 
  
ggplot(niche_data_not_scaled, aes(x = value, fill = name)) +
    geom_density(alpha = 0.5, cut = 0) +
    facet_wrap(~name, scales = "free", nrow = 1) +
    labs(title = "Pre-extinction",
         x = NULL,
         y = "Occurrence Density") +
    theme_bw() +
    theme (legend.position = "none",
           plot.title = element_text(hjust = 0.5, face = "bold")
    ) 

combined_niche <- bind_rows(
  niche_data %>% mutate(Period = "Pre-extinction"),
  niche_data_post %>% mutate(Period = "Post-extinction")
) %>%
  mutate(Period = factor(Period, levels = c("Pre-extinction", "Post-extinction")))

make_var_plot <- function(var_name, color_hex) {
  df_sub <- combined_niche %>% filter(Variable == var_name)
  
  display_names <- c("DOM" = "DOM", "H2S" = expression(H[2]*S), 
                     "O2" = expression(O[2]), "pH" = "pH", 
                     "SO4" = expression(SO[4]), "temp" = "Temperature (°C)")
  
  p <- ggplot(df_sub, aes(x = Value, fill = Period)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, alpha = 0.6, position = "identity", color = "white", size = 0.1) +
    geom_density(size = 0.8, alpha = 0.2) +
    facet_grid(Period ~ .) + 
    scale_fill_manual(values = c("Pre-extinction" = color_hex, "Post-extinction" = color_hex)) +
    labs(title = display_names[[var_name]], x = NULL, y = NULL) +
    theme_bw(base_size = 14) + 
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text = element_text(size = 10, color = "black"),
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(size = 11, face = "bold") 
    )
  
  if(var_name != "temp") {
    p <- p + theme(strip.text.y = element_blank(), 
                   strip.background = element_blank())
  }
  
  return(p)
}
p_dom  <- make_var_plot("DOM", "#F8766D")
p_h2s  <- make_var_plot("H2S", "#B79F00")
p_o2   <- make_var_plot("O2",  "#00BA38")
p_ph   <- make_var_plot("pH",  "#00BFC4")
p_so4  <- make_var_plot("SO4", "#619CFF")
p_temp <- make_var_plot("temp", "#F564E3")

final_layout <- (p_dom | p_h2s | p_o2 | p_ph | p_so4 | p_temp) + 
  theme(
    plot.subtitle = element_text(size = 12, hjust = 0.5, face = "italic"),
    plot.margin = margin(10, 10, 10, 10)
  )

final_layout

niche_summary <- combined_niche %>%
  group_by(Variable, Period) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE),
    Breadth_IQR = IQR(Value, na.rm = TRUE), 
    Sample_Size = n()
  ) %>%
  ungroup()

################
genus_counts <- epme_data %>%
  group_by(genus) %>%
  summarise(n_points = n(), .groups = "drop")

max_val <- max(genus_counts$n_points)
breaks_seq <- seq(0, max_val + 10, by = 10)

genera_stats <- genus_counts %>%
  mutate(bin = cut(n_points, 
                   breaks = breaks_seq, 
                   include.lowest = TRUE,
                   dig.lab = 5)) %>% # Ensures large numbers aren't in scientific notation
  group_by(bin) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(percentage = (count / sum(count)) * 100) %>%
  ungroup()

######################## Motility density
motility_list <- unique(niche_data$motility)
plot_list <- list()

for (i in 1:length(motility_list)) {
  p <- motility_list[i]
  temp_data <- niche_data %>% filter(motility == p)
  
  p_plot <- ggplot(temp_data, aes(x = Value, fill = Variable)) +
    geom_histogram(aes(y = after_stat(density)), 
                   bins = 30, alpha = 0.8, position = "identity") +
    geom_density(alpha = 0.4, colour = "black") +
    # CHANGE 1: Use ncol = 1 and strip.position = "left" to stack variables vertically
    facet_wrap(~Variable, scales = "free", ncol = 1, strip.position = "left") +
    scale_fill_manual(values = my_colors) +
    labs(title = p, x = NULL, y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside", # Puts labels outside the axis
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 8)
    )
  
  if (i == 1) {
    p_plot <- p_plot + theme(
      strip.text.y = element_text(angle = 0, face = "bold", size = 10)
    )
  } else {
    p_plot <- p_plot + theme(
      strip.text.y = element_blank()
    )
  }
  
  plot_list[[i]] <- p_plot
}

final_stitch <- plot_list[[1]] | plot_list[[2]] | plot_list[[3]]

final_stitch

############################
auc_over_7_pre <- ENM_results_pre %>% filter(test_AUC >= 0.7 & test_AUC < 1)
auc_over_7_post <- ENM_results_post %>% filter(test_AUC >= 0.7 & test_AUC < 1)

auc_over_7_both <- intersect(auc_over_7_pre$genus, auc_over_7_post$genus)

niche_summary <- combined_importance %>%
  group_by(Variable, Epoch) %>%
  summarise(
    Mean = mean(Importance, na.rm = TRUE),
    Median = median(Importance, na.rm = TRUE),
    SD = sd(Importance, na.rm = TRUE),
    Min = min(Importance, na.rm = TRUE),
    Max = max(Importance, na.rm = TRUE),
    # Niche Breadth (Interquartile Range) is great for reporting
    Breadth_IQR = IQR(Importance, na.rm = TRUE), 
    Sample_Size = n()
  ) %>%
  ungroup()