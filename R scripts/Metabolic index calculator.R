install.packages("ncdf4")
library(ncdf4)

epme_data <- read.csv("C:/Users/am14g22/OneDrive - University of Southampton/Paleodistribution paper/epme_data.csv")

# pre and post
epme_data_pre <- epme_data %>% filter(stg <= 51)
epme_data_post <- epme_data %>% filter(stg > 51)

# Load climate model experiment 1 and 3
exp1 <- "C:/Users/am14g22/OneDrive - University of Southampton/Paleodistribution paper/HulseNatGeo2020/HulseNatGeo2020/01_endPermian_background_02.5CO2_1.00PO4_Tmp_98sink_sulfE6"
exp3 <- "C:/Users/am14g22/OneDrive - University of Southampton/Paleodistribution paper/HulseNatGeo2020/HulseNatGeo2020/03_MainExtinction_20.0CO2_2.00PO4_Tmp_98sink_sulfE6"

# Load environmental data for experiment 1 (pre extinciton)
lvl2_O2_1 <- cGENIE.data(var = "ocn_O2", experiment = exp1, depth.level = 2)
lvl2_temp_1 <- cGENIE.data(var = "ocn_temp", experiment = exp1, depth.level = 2)
lvl2_sal_1 <- cGENIE.data(var = "ocn_sal", experiment = exp1, depth.level = 2)

# Load environmental data for experiment 3 (post extinction)
lvl2_O2_3 <- cGENIE.data(var = "ocn_O2", experiment = exp3, depth.level = 2)
lvl2_temp_3 <- cGENIE.data(var = "ocn_temp", experiment = exp3, depth.level = 2)
lvl2_sal_3 <- cGENIE.data(var = "ocn_sal", experiment = exp3, depth.level = 2)

# Define datasets with important variables (pre and post)
pre_env_dataset_metabolic_index <- list(
  O2 = lvl2_O2_1, temp = lvl2_temp_1, sal = lvl2_sal_1)

post_env_datasets_metabolic_index <- list(
  O2 = lvl2_O2_3, temp = lvl2_temp_3, sal = lvl2_sal_3)

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
  resolution = c(1, 1),
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
for (var_name in names(pre_env_dataset_metabolic_index)) {
  message(paste("Extracting:", var_name))
  
  extractor <- nearest_ocean_extractor(pre_env_dataset_metabolic_index[[var_name]], max_dist = tolerate_dist)
  epme_data_pre[[var_name]] <- sapply(1:nrow(epme_data_pre), function(i) {
    extractor(epme_data_pre$paleolng[i], epme_data_pre$paleolat[i])
  })
}

# Remove NAs from points 
epme_data_pre <- epme_data_pre %>% 
  filter(!is.na(temp))

# Generate rasters
raster_list <- lapply(pre_env_dataset_metabolic_index, function(df) {
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
for (var_name in names(post_env_datasets_metabolic_index)) {
  message(paste("Extracting:", var_name))
  
  extractor <- nearest_ocean_extractor(post_env_datasets_metabolic_index[[var_name]], max_dist = tolerate_dist)
  epme_data_post[[var_name]] <- sapply(1:nrow(epme_data_post), function(i) {
    extractor(epme_data_post$paleolng[i], epme_data_post$paleolat[i])
  })
}

# Remove NAs from points 
epme_data_post <- epme_data_post %>% 
  filter(!is.na(temp))

# Generate rasters
raster_list <- lapply(post_env_datasets_metabolic_index, function(df) {
  cGENIE_df_to_raster(
    df = df,
    value_col = "var",
    template_rast = template_rast
  )
})

# Convert the list to a multi-layer SpatRaster object
Env_rasters_all_post <- rast(raster_list)
Env_rasters_all_post <- mask(Env_rasters_all_post, sum(Env_rasters_all_post))

#Define functions needed to calculate MI

#press
press <- function(d, lat)
{
  press <- c()
  for (i in 1:length(d))
  {
    press <- c(press, gauge_p(d=d[[i]], lat=lat[[i]]))
  }
  return(press)
}

#p_corr_pp
p_corr_pp <- function(pp, p, t, species="O2")
{
  nu <- switch(species,
               O2 =mean(31.9, 31.9, 32.2, 32.3, 32, 31.7),      # cm^3/mol (molar volume of O2 here)
               CO2=mean(35.3, 34.4),
               N2 =mean(33.2, 33.2, 33.5),
               Ar =mean(32.3, 32.1),
               He =mean(29.7, 29.7))
  R  <- 8.314472  # J/(mol*K) = Nm/(mol*K) = kg*m^2/s^2/(mol*K)
  P1 <- 0         # atm = 1.01325 bar = 1.01325*10^5 N/m^2
  T  <- t+273.15  # K
  p  <- p/1.01325 # atm = 1.01325 bar = 1.01325*10^5 N/m^2
  
  return(exp(-(nu*(P1-p)/(R*T*10)))*pp)
}

#pO2_ multiply o2 by 1,000,000 to get umol from mol
pO2_ <- function(O2, t, S, p)
{
  o2 <- o2 * 1000000
  thet    <- theta(S=S, t=t, p=p)      # deg C
  O2_sat  <- O2sat_new(S=S, t=thet)    # umolO2/kg-soln
  xO2air  <- 0.20946                   # mole fraction: dimensionless
  
  pO2  <- (O2/O2_sat)*(xO2air*(1-sat_vap_press(S=S, t=thet))*1e6)   # at surface: total pressure = 1 atm 
  attr(pO2, "unit") <- "uatm"
  
  return(pO2)
}

#theta
theta <- function(S, t, p)
{
  PR <- 0                          # reference pressure
  H  <- PR - p                     # pressure differential
  
  tpar  <- c(0.5, 0.29289322,   1.707106781, 0)
  qpar1 <- c(1,   0.58578644,   3.414213562, 0)
  qpar2 <- c(0,   0.121320344, -4.121320344, 1)
  
  Q <- 0; XK <- 0 
  for (i in 1:4)
  {
    XK <- H*delta_adiabat(S, t, p)
    
    t  <- t + tpar[[i]]*(XK - Q)
    Q  <- qpar1[[i]]*XK + qpar2[[i]]*Q
    
    if (!(i==2)) {p  <- p + 0.5*H}
  }
  theta = t + (XK - 2*Q)/6
  
  attr(theta, "unit") <- "degC"
  return(theta)
}

#delta_adiabat
delta_adiabat <- function(S, t, p)
{
  DS <- S - 35    
  p  <- p*10                            # convert bar to dbar
  
  delta_adiabat <- ((((-2.1687e-16*t + 1.8676e-14)*t - 4.6206e-13)*p + ((2.7759e-12*t - 1.1351e-10)*DS
                                                                        +((-5.4481e-14*t + 8.733e-12)*t - 6.7795e-10)*t + 1.8741e-8))*p
                    +(-4.2393e-8*t + 1.8932e-6)*DS+((6.6228e-10*t - 6.836e-8)*t + 8.5258e-6)*t + 3.5803e-5)
  
  delta_adiabat <- delta_adiabat*10     # convert /dbar to /bar
  
  attr(delta_adiabat, "unit") <- "degC/bar"
  return(delta_adiabat)
}

#O2sat_new
O2sat_new <- function(S, t)
{
  #'CALCULATE OXYGEN CONCENTRATION AT SATURATION
  #'Source:  Oxygen solubility in seawater: Better fitting equations.
  #'          Garcia & Gordon (1992) L&O 37: 1307-1312.
  #'Input:       S = Salinity (pss-78)
  #'             T = Temp (deg C)
  #'Output:      Oxygen saturation at one atmosphere (umol/kg).
  
  #'DEFINE CONSTANTS, ETC FOR SATURATION CALCULATION
  
  #' The constants used are for units of umol O2/kg.
  T <- t
  
  A0 = 5.80871
  A1 = 3.20291
  A2 = 4.17887
  A3 = 5.10006
  A4 = -0.0986643
  A5 = 3.80369
  
  B0 = -0.00701577
  B1 = -0.00770028
  B2 = -0.0113864
  B3 = -0.00951519
  
  C0 = -0.000000275915
  
  #'CALCULATE Ts FROM T (deg C):
  
  TS = log((298.15 - T) / (273.15 + T))
  
  #'CALCULATE O2 SATURATION in umol/kg:
  
  A = ((((A5 * TS + A4) * TS + A3) * TS + A2) * TS + A1) * TS + A0
  B = ((B3 * TS + B2) * TS + B1) * TS + B0
  O2sat_new = exp(A + S * (B + S * C0))
}

#sat_vap_press
sat_vap_press <- function(S, t)
{
  T = t + 273.15
  return(exp(24.4543 - (6745.09/T) - 4.8489*log(T/100) - 0.000544*S))
}

#pO2
pO2 <- function(O2, S, t, d, lat) 
  
{ 
  
  p   <- press(d=d, lat=lat) 
  
  pO2 <- (p_corr_pp(pp=pO2_(O2=O2, t=t, S=S, p=p), p=p, t=t))/1000 
  
  attr(pO2, "unit") <- "matm" 
  
  return(pO2) 
  
}

#Add depth column for ease 
epme_data_pre$depth <- 127.552
epme_data_post$depth <- 127.552

# work out pO2 for pre extinction occurences
epme_data_pre$pO2 <- mapply(
  FUN  = pO2,
  O2   = epme_data_pre$O2,
  S    = epme_data_pre$sal,
  t    = epme_data_pre$temp,
  d    = epme_data_pre$depth,
  lat  = epme_data_pre$paleolat
)

# work out pO2 for post extinction occurences
epme_data_post$pO2 <- mapply(
  FUN  = pO2,
  O2   = epme_data_post$O2,
  S    = epme_data_post$sal,
  t    = epme_data_post$temp,
  d    = epme_data_post$depth,
  lat  = epme_data_post$paleolat
)

Ao_and_Eo_averages <- Ao_and_Eo %>%
  group_by(`phylum`) %>%
  summarise(
    mean_Ao = mean(Ao, na.rm = TRUE),
    mean_Eo = mean(Eo, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

epme_data_pre <- epme_data_pre %>%
  left_join(
    Ao_and_Eo_averages %>% select("phylum", mean_Ao, mean_Eo),
    by = "phylum"
  )

metabolic_index <- function(Ao, pO2, Eo, t, t_ref)
{
  pO2_atm <- pO2 / 1000
  kB = 8.617333262e-5
  t_k <- t + 273.15
  t_ref_k <- t_ref + 273.15
  
  phi<- Ao * pO2_atm / exp((-Eo / kB) * (1/t_k - 1/t_ref_k))
  
  return(phi)
    
}

epme_data_pre$phi <- metabolic_index(
  Ao = epme_data_pre$mean_Ao,
  pO2 = epme_data_pre$pO2,
  Eo = epme_data_pre$mean_Eo,
  t = epme_data_pre$temp,
  t_ref = 15)
 

# Find the 3D output file (usually biogem/fields_biogem_3d.nc)
list.files(exp1, pattern = "3d.*\\.nc$", recursive = TRUE, full.names = TRUE)

nc_file <- file.path(exp1, "biogem", "fields_biogem_3d.nc")   # adjust to what the line above prints
nc <- nc_open(nc_file)
zt <- ncvar_get(nc, "zt")        # mid-point depth of each layer, m
nc_close(nc)

zt[2]
