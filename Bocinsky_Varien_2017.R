## This is the script for the analyses in Varien and Bocinsky, under review:
## R. Kyle Bocinsky and Mark D. Varien. 2017. Calibrating Maize 
## Paleoproduction Models using Experimental Data. Ethnobiology. 37(2).

## The analysis consists of three components. First, VEP-style paleoproduction
## estimates are generated for the locations of the Pueblo Farming Project gardens.
## Then, PFP data are analysed to estimate standard productivity across gardens.
## Finally, these data are brought together in a series of figures to be used in
## the manuscript.

## This script makes heavy use of piping from the magrittr package 
## and the data manipulation functions in the dplyr package. See 
## https://cran.r-project.org/web/packages/magrittr/vignettes/magrittr.html
## and
## https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html

## Throughout the VEP productivity reconstruction, references are made to 
## the canonical VEP reconstruction reference: Kohler 2012

## Kohler, Timothy A. 2012. Modeling agricultural productivity and farming effort.
## In Kohler, T. A. and Varien, M. D., editors, Emergence and Collapse of Early 
## Villages: Models of Central Mesa Verde Archaeology, chapter 6, pages 85–112.
## University of California Press, Berkeley, California.

## Author: R. Kyle Bocinsky
## Date: 1/3/2017

## Load requisite packages
library(FedData)
FedData::pkg_test("magrittr")
FedData::pkg_test("Hmisc")
FedData::pkg_test("lubridate")
FedData::pkg_test("readr")
FedData::pkg_test("zoo")
FedData::pkg_test("stringr")
FedData::pkg_test("sp")
FedData::pkg_test("raster")
FedData::pkg_test("rasterVis")
FedData::pkg_test("ggplot2")
FedData::pkg_test("gridExtra")
FedData::pkg_test("raster")
FedData::pkg_test("sp")
FedData::pkg_test("RColorBrewer")
FedData::pkg_test("rgeos")
FedData::pkg_test("broom")
FedData::pkg_test("maptools")
FedData::pkg_test("dplyr")
FedData::pkg_test("doParallel")
FedData::pkg_test("geojsonio")

# Load other useful functions
for(f in list.files("./src", pattern = ".R", full.names = T)){
  source(f)
}

# Suppress scientific notation
options(scipen=999)

# Perform a "clean" run?
clean = T
if(clean){
  unlink("./OUTPUT", recursive = T, force = T)
  unlink("./Bocinsky_Varien_2017.zip", force = T)
}

# Create directories for output
dir.create("./OUTPUT/DATA", recursive = T, showWarnings = F)
dir.create("./OUTPUT/FIGURES", recursive = T, showWarnings = F)
dir.create("./OUTPUT/TABLES", recursive = T, showWarnings = F)

# Growing seasons of the PFP
seasons <- 2009:2015


##### BEGIN VEP PRODUCTIVITY ESTIMATION #####

### SOILS & PRODUCTIVITY ###
# Import VEPI Paleoproductivity
vepi.paleoprod <- raster::brick("./DATA/vepi_paleoprod.tif")

# Load a polygon of the Indian Camp Ranch + Crow Canyon study area
ccac <- rgdal::readOGR(dsn = "./DATA/ccac.geojson", "OGRGeoJSON", verbose = FALSE)

# Load the VEP I soils data for Indian Camp Ranch / CCAC
vepi.soils.ccac <- rgdal::readOGR(dsn = "./DATA/vepi_soils_ccac.geojson", "OGRGeoJSON", verbose = FALSE)

### PDSI ###
# Kohler 2012: Steps 2--3 (pp. 90--92)
# When to start the burn-in for PDSI
burnin_start <- 1991

# What years to calculate PDSI norms over
pdsi_norm_years <- 1981:2010

# Use FedData to download the weather data from the Cortez GHCN
# weather station, minimum and maximum temperature and precipitation
# Kohler 2012: Step 1 (pp. 86--90)
cortez_weather <- c(FedData::get_ghcn_daily_station(ID="USC00051886", 
                                                    elements = c("TMIN","TMAX"), 
                                                    standardize = T, 
                                                    raw.dir = "./OUTPUT/DATA/GHCN"),
                    FedData::get_ghcn_daily_station(ID="USC00051886", 
                                                    elements = c("PRCP"), 
                                                    raw.dir = "./OUTPUT/DATA/GHCN"))

# Calculate monthly weather data from the daily GHCN data
cortez_weather_monthly <- cortez_weather %>%
  FedData::station_to_data_frame() %>% # convert station data matrix to a data frame
  dplyr::as_data_frame() %>% # convert to a tibble
  dplyr::filter(lubridate::year(DATE) %in% burnin_start:max(seasons), # Get only rows with dates after the burn-in date, and during the PFP seasons
                DATE < lubridate::floor_date(Sys.Date(), unit = "month")) %>%  # and only rows with dates prior to the current month
  dplyr::mutate(DATE = lubridate::ymd(DATE), # convert dates to Date objects
                TMIN = zoo::na.approx(TMIN/10, na.rm = F), # Fill NA values in TMIN by linear interpolation
                TMAX = zoo::na.approx(TMAX/10, na.rm = F), # Fill NA values in TMAX by linear interpolation
                TMAX = ((TMAX)*1.8 + 32), # Convert TMAX to Fahrenheit
                TMIN = ((TMIN)*1.8 + 32), # Convert TMIN to Fahrenheit
                TAVG = (TMAX+TMIN)/2, # Calculate the average temperature
                PRCP = PRCP*0.00393701) %>% # Convert precipitation to hundreths of an inch
  dplyr::rename(TAVG_F = TAVG, PRCP_IN = PRCP) %>% # Change names to make them more explanatory
  dplyr::mutate(MONTH = lubridate::month(DATE), YEAR = lubridate::year(DATE)) %>% # Extract months and years from dates
  dplyr::select(YEAR,MONTH,TAVG_F,PRCP_IN) %>% # Keep only these variables
  dplyr::group_by(YEAR,MONTH) %>% # Group the data_frame by year and month
  dplyr::summarise(TAVG_F = mean(TAVG_F), PRCP_IN = sum(PRCP_IN, na.rm = T)) # Calculate TAVG averages and net precip per month

# Calculate monthly norms over 1981 -- 2010
cortez_weather_monthly_norms <- cortez_weather %>%
  station_to_data_frame() %>% # convert station data matrix to a data frame
  dplyr::as_data_frame() %>% # convert to a tibble
  dplyr::filter(lubridate::year(DATE) %in% pdsi_norm_years) %>% # Get only rows with years in the pdsi_norm_years
  dplyr::mutate(DATE = lubridate::ymd(DATE), # convert dates to Date objects
                TMIN = zoo::na.approx(TMIN/10, na.rm = F), # Fill NA values in TMIN by linear interpolation
                TMAX = zoo::na.approx(TMAX/10, na.rm = F), # Fill NA values in TMAX by linear interpolation
                TMAX = ((TMAX)*1.8 + 32), # Convert TMAX to Fahrenheit
                TMIN = ((TMIN)*1.8 + 32), # Convert TMIN to Fahrenheit
                TAVG = (TMAX+TMIN)/2, # Calculate the average temperature
                PRCP = PRCP*0.00393701) %>% # Convert precipitation to hundreths of an inch
  dplyr::rename(TAVG_F = TAVG, PRCP_IN = PRCP) %>% # Change names to make them more explanatory
  dplyr::mutate(MONTH = lubridate::month(DATE), YEAR = lubridate::year(DATE)) %>%  # Extract months and years from dates
  dplyr::select(YEAR,MONTH,TAVG_F,PRCP_IN) %>% # Keep only these variables
  dplyr::group_by(YEAR,MONTH) %>% # Group the data_frame by year and month
  dplyr::summarise(TAVG_F = mean(TAVG_F), PRCP_IN = sum(PRCP_IN, na.rm = T)) %>% # Calculate TAVG averages and net precip per month
  dplyr::ungroup() %>% #ungroup
  dplyr::select(MONTH,TAVG_F,PRCP_IN) %>% # Keep only these variables
  dplyr::group_by(MONTH) %>% # Group the data_frame by month
  dplyr::summarise(TAVG_F = mean(TAVG_F), PRCP_IN = mean(PRCP_IN)) # Calculate averages per month

# Fill missing months with norms by first calculating all months that should exist,
# and then finding the difference with those that do exist.
# April and May of 1997 were missing.
cortez_weather_monthly <- cortez_weather_monthly %>% 
  dplyr::bind_rows(
    anti_join(
      cortez_weather_monthly_norms %>%
        dplyr::full_join(
          expand.grid(YEAR = burnin_start:max(seasons), MONTH = 1:12) %>% # Get full combination of years and months
            filter(
              lubridate::ymd(paste0(sprintf("%04d", YEAR),sprintf("%02d", MONTH),"01")) < lubridate::floor_date(Sys.Date(), unit = "month")), # filter only months up to the present month
          by = c("MONTH")), # Join with cortez_weather_monthly_norms by month
      cortez_weather_monthly, # Anti-join with cortez_weather_monthly by month and year
      by=c("YEAR","MONTH")
      )
    ) %>%
  dplyr::arrange(YEAR,MONTH) %>% # Sort by year and month
  dplyr::full_join(expand.grid(YEAR = burnin_start:max(seasons), MONTH = 1:12), by=c("YEAR","MONTH")) # Add NAs at end if ending in current year

# Get the VEP I latitude for calculating PDSI
# Following Kohler 2012:92 we use
# the latitude at the center of the study area
vepi_latitude <- vepi.paleoprod %>%
  FedData::polygon_from_extent() %>%
  sp::spTransform(sp::CRS("+proj=longlat +ellps=WGS84")) %>%
  sp::bbox() %>%
  rowMeans() %>%
  magrittr::extract("y")

## Produce PDSI estimates for each soil
# Download the scpdsi C code from GitHub
FedData::download_data("https://raw.githubusercontent.com/crowcanyon/vep_paleoprod/master/src/scpdsi_r.cpp",
                       destdir = "./src/")

# Load the PDSI c program
Rcpp::sourceCpp("./src/scpdsi_r.cpp")

# Run the PDSI program for each soil
# This operates in parallel; it will register the number of cores appropriate to your machine
registerDoParallel()
vepi.soils.ccac.pdsi <- foreach::foreach(soil = 1:nrow(vepi.soils.ccac)) %dopar% {
  
  out <- scpdsi(monthly_T = cortez_weather_monthly$TAVG_F,
         monthly_P = cortez_weather_monthly$PRCP_IN,
         mon_T_normal = cortez_weather_monthly_norms$TAVG_F,
         awc = vepi.soils.ccac@data[soil,"AWC_Lower_median"],
         lat = vepi_latitude,
         start_year = burnin_start,
         end_year = max(seasons))
  expand.grid(YEAR = burnin_start:max(seasons), MONTH = 1:12) %>%
    tibble::as_tibble() %>%
    dplyr::arrange(YEAR,MONTH) %>%
    dplyr::mutate(PDSI = out) %>%
    dplyr::filter(YEAR %in% seasons, MONTH == 6) %>%
    dplyr::select(YEAR,PDSI) 
  
}
stopImplicitCluster()

### RECONSTRUCTION ###
## There are 5 components to the reconstruction back in time: 
# The measured historical maize yields collected by Burns
# The calculated June PDSI values for "Bean Soils"
# The MV Douglas Fir tree-ring chronology (proxying PDSI)
# The 1st Principle Component of the combined San Francisco Peaks and Almagre chronologies (proxying September mean temperature)

## Here, to perform a modern "reconstruction" we use the calculated PDSI directly.
## Read in the calibration data used in VEP I
## These include Burns' estimates of maize yields in bushels/acre, 
## and the PDSI estimates for bean soils
VEPI.calibration.data <- readr::read_csv("./DATA/VEPI_calibration_data.csv")

# Kohler 2012: Step 5 (pp. 92--95)
# Calculate the linear models of the Burns data on PDSI, controlling for year
# YEAR is included as an explicit predictor to correct for a "technology trend" in 1931--1960
VEPI_LM <- lm(MAIZE ~ scale(YEAR) + PDSI, data=VEPI.calibration.data)
pdsi.coef <- coefficients(VEPI_LM)['PDSI']
intercept <- coefficients(VEPI_LM)['(Intercept)']

## Estimate potential maize production for each year, for each soil
# Perform the prediction
# Retrodiction equation (per cell): prediction ~ pdsi.coef*modern.PDSI + intercept
predictions <- lapply(vepi.soils.ccac.pdsi,function(x){
  ((pdsi.coef * x$PDSI) + intercept) * 62.77 # Convert to kg/ha
})
predictions <- do.call(rbind,predictions)

## Re-weight production for each soil 
# Read in the NRCS normal year dry-weight soil productivity (NPP) & Bean Soils
NPP.bean.mean <- 1093 # Mean Bean Soil NPP for all VEP I soils, in lbs/ac
NPP.reweight <- vepi.soils.ccac$NYProd_lb_ac/NPP.bean.mean
predictions <- sweep(predictions, MARGIN=1, NPP.reweight,`*`)

## We depart a bit from the methods described in Kohler 2012. Before, the retrodicted values
## would be multiplied by the NPP renormed values, then loaded into the simulation, where they 
## would be renormed to prehispanic maize varieties and further reduced by a cold correction 
## and a hand-planting reduction (steps 9, 10, and 11 in Kohler 2012).
## We perform all of those steps "up front" here.

## RENORM MAIZE PRODUCTION FOR PREHISPANIC VARIETIES AND CULTIVATION PRACTICES
# Kohler 2012: Step 9 (pp. 100--103)
renorm.factor <- 0.68
predictions <- predictions * renorm.factor

## MAKE HAND-PLANTING ADJUSTMENT
# Kohler 2012: Step 10 (pp. 103)
# Read in hand-planting factor from the soils data
predictions <- sweep(predictions,MARGIN=1,vepi.soils.ccac$SCM_RED,`*`)
colnames(predictions) <- seasons

## MAKE COLD CORRECTION
# Kohler 2012: Step 11 (pp. 103--106)
# We don't make the cold correction here, because CCAC PFP gardens are below 7,000 feet.

# Write yield info to soils geojson file
vepi.soils.ccac@data <- cbind(vepi.soils.ccac@data,predictions %>% tibble::as_tibble())
geojsonio::geojson_write(vepi.soils.ccac,file = "./OUTPUT/DATA/vepi.soils.ccac_yields.geojson")

# Rasterize data based on VEPI 200m simulation grid
# A function to rasterize each variable in a Spatial* object, and write the outputs as GeoTIFFs
# This operates in parallel; it will register the number of cores appropriate to your machine
rasterize_each <- function(x, y, force.redo = F){
  # Transform the polygons into the coordinate reference system of the study area
  x %<>%
    sp::spTransform(sp::CRS(raster::projection(y)))
  
  registerDoParallel()
  out <- foreach::foreach(var = names(x), .combine = raster::stack) %dopar% {
    if(file.exists(paste0("./OUTPUT/DATA/vepi_soils_ccac_raster_",var,".tif")) & !force.redo)
      return(raster::raster(paste0("./OUTPUT/DATA/vepi_soils_ccac_raster_",var,".tif")))
    
    out <- x %>%
      raster::rasterize(y = y,
                        field = var,
                        fun = mean) %>%
      raster::crop(y = x, snap = "out") %T>%
      writeRaster(paste0("./OUTPUT/DATA/vepi_soils_ccac_raster_",var,".tif"),
                  datatype = "FLT4S",
                  options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                  overwrite = T,
                  setStatistics = FALSE)
    return(out)
  }
  stopImplicitCluster()
  return(out)
}

# Rasterize and write each variable
vepi.soils.ccac.raster <- vepi.soils.ccac[,c("CLUSTER","VPSCode","AWC_Lower_median","NYProd_lb_ac","SCM_RED",as.character(seasons))] %>%
  rasterize_each(y = vepi.paleoprod)
names(vepi.soils.ccac.raster) <- c("CLUSTER","VPSCode","AWC_Lower_median","NYProd_lb_ac","SCM_RED",as.character(seasons))
##### END VEP PRODUCTIVITY ESTIMATION #####

##### BEGIN PFP PRODUCTIVITY ESTIMATION #####
gardens <- readr::read_csv("./DATA/gardens.csv")

# Read in the garden locations and metadata
garden_locations <- rgdal::readOGR(dsn = "./DATA/gardens.geojson", "OGRGeoJSON", verbose = FALSE)

# Read in data recorded about each ear harvested
ears <- readr::read_csv("./DATA/ears.csv")

# Estimate Kernel and Cob weight for ears withheld whole from analysis (e.g., POG, 2009, Clump 24)
# This is done by calculating the average ratio of kernel to ear weight, and then extrapolating to
# cob and kernel weight for samples for which we only have ear weight.
ears %<>%
  dplyr::group_by(Season,Garden,Condition) %>% # Group by season, garden, and condition
  dplyr::summarise(Kernel_r = mean(`Kernel weight`/`Ear weight`, na.rm = T)) %>% # Calculate kernel to ear ratio
  dplyr::mutate(Cob_r = 1 - Kernel_r) %>% # define cob to ear ratio as 1 minus the kernel to ear ratio
  dplyr::full_join(ears, by = c("Season","Garden","Condition")) %>% # join results back to ears tibble
  dplyr::mutate(`Kernel weight` = ifelse(is.na(`Kernel weight`),`Ear weight` * Kernel_r,`Kernel weight`), # Fill missing values
                `Cob weight` = ifelse(is.na(`Cob weight`),`Ear weight` * Cob_r,`Cob weight`)) %>%
  dplyr::select(-Kernel_r, -Cob_r) %>% # drop the ratio columns
  dplyr::ungroup() %>%
  dplyr::left_join(y = (gardens %>% dplyr::select(Season,Garden,Variety)), by = c("Season","Garden")) %>% # join with the gardens data
  dplyr::mutate(Variety = as.factor(Variety)) %>% # turn "variety" into a categorical variable
  dplyr::arrange(Season, Garden, Clump) %>% # sort by these variables
  dplyr::select(Season, Garden, Variety, Clump, Condition, Rows, `Ear weight`, `Cob weight`, `Kernel weight`) # reorder columns

# Create a grid of expected clumps, if all clumps were weighed (they weren't)
expected.clumps <- gardens %>%
  dplyr::select(Season, Garden, Clumps) %>%
  split(list(gardens$Season, gardens$Garden)) %>%
  lapply(function(x){
    expand.grid(Season = x$Season, Garden = x$Garden, Clump = 1:x$Clumps)
  }) %>%
  bind_rows()

# Estimate kernel yields by simulating distributions across clump yields
# Here, we calculate three yields:
# the first is just the "raw" yield of kernel weight per garden area;
# in the second, we only calculate over planted area as defined by an estimate of "spacing" between clumps
# in the third, we standardize the density of planting to 2m plant spacing (4 sq m per clump; Beaglehole 1937:40; Bellorado 2007:96; Bradfield 1971:5; Dominguez and Kolm 2003, 2005)
yields <- ears %>%
  dplyr::select(Season, Garden, Clump, `Kernel weight`) %>% #select these columns
  dplyr::group_by(Season, Garden, Clump) %>%# calculations are by season and garden
  dplyr::summarise(`Net kernel weight` = sum(`Kernel weight`)) %>%
  dplyr::full_join(expected.clumps, by = c("Season","Garden","Clump")) %>%
  dplyr::left_join(y = gardens, by = c("Season","Garden")) %>% # join back to gardens
  dplyr::ungroup() %>%
  dplyr::mutate(`Net kernel weight` = ifelse(is.na(`Net kernel weight`),0,`Net kernel weight`)) %>% #recode missing values to zeros; those clumps didn't produce
  dplyr::select(Season:Clumps, Clumps, Spacing) %>% #select these columns
  dplyr::mutate(#`Yield by clump area` = (Clumps * `Net kernel weight`/1000)/((Spacing ^ 2) * Clumps * 0.0001), # yield by actual clump area
                `PFP experimental yield` = (Clumps * `Net kernel weight`/1000)/((2 ^ 2) * Clumps * 0.0001) # yield by clump area with 2m spacing
                ) %>% 
  dplyr::mutate(Variety = as.factor(Variety),
                Garden = as.factor(Garden),
                Season = as.factor(Season)) %>%
  dplyr::arrange(Season, Garden, Clump) %>% # sort by these variables
  dplyr::select(Season,
                Garden,
                Variety,
                Clumps,
                Spacing,
                Clump,
                `Net kernel weight`,
                # `Yield by clump area`,
                `PFP experimental yield`) %>% # reorder columns
  dplyr::group_by(Season, Garden) %>% # calculations are by season and garden 
  dplyr::rename(`Spacing (m)` = Spacing,
                `Net kernel weight (g)` = `Net kernel weight`,
                `PFP experimental yield (kg/ha)` = `PFP experimental yield`)

readr::write_csv(yields,"./OUTPUT/DATA/yields.csv")

##### END PFP PRODUCTIVITY ESTIMATION #####

##### BEGIN JOINT ANALYSIS #####
## Extract the VEP yield reconstructions under the PFP gardens
# Read in the VEP productivity reconstruction from above
soils_VEPI_yields <- rgdal::readOGR(dsn = "./OUTPUT/DATA/vepi.soils.ccac_yields.geojson", "OGRGeoJSON", verbose = FALSE)
soils_VEPI_yields$muname <- gsub(" MLRA 36","",soils_VEPI_yields$muname)

PFP_VEP_yields <- garden_locations %>%
  rgeos::gCentroid(byid=T) %>% # Get the centroid under each garden
  sp::spTransform(sp::CRS(raster::projection(vepi.soils.ccac.raster))) %>% # transform into the same projection
  raster::extract(x = vepi.soils.ccac.raster) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(Garden = garden_locations$Abbreviation) %>%
  tidyr::gather(Season, Yield, X2009:X2015) %>%
  dplyr::as_data_frame() %>%
  dplyr::mutate(Season = gsub("X","",Season) %>% as.factor()) %>%
  dplyr::select(Garden,Season,Yield) %>%
  dplyr::arrange(Garden, Season)

##### garden_soils_map.pdf #####
## A map showing the soils on the CCAC campus and the locations of the PFP gardens
## Map of NRCS soil complexes (mapunits) on the Crow Canyon campus, 
## and the locations of experimental gardens reported here. 
## CDG: check dam garden; KUG: Karen’s upper garden; 
## PLC: Pueblo Learning Center garden; POG: Paul’s old garden.
# Prepare the data for plotting in ggplot2
soils_VEPI_yields %<>% sp::spTransform(sp::CRS(raster::projection(vepi.soils.ccac.raster)))
soils_VEPI_yields@data$id = rownames(soils_VEPI_yields@data)
soils_VEPI_yields.points = ggplot2::fortify(soils_VEPI_yields, region="id") %>% as_data_frame()
soils_VEPI_yields.df = dplyr::full_join(soils_VEPI_yields.points, soils_VEPI_yields@data, by="id")
garden_locations %<>% sp::spTransform(sp::CRS(raster::projection(vepi.soils.ccac.raster)))
garden_locations@data$id = rownames(garden_locations@data)
garden_locations.points = ggplot2::fortify(garden_locations, region="id") %>% as_data_frame()
garden_locations.df = dplyr::full_join(garden_locations.points, garden_locations@data, by="id")

vep_soil_polygons <- ggplot(soils_VEPI_yields.df) +
  aes(long,lat,group=group,fill=factor(VPSCode, levels = sort(unique(VPSCode)))) +
  geom_polygon() + 
  ggtitle("Soil mapunit polygons") +
  geom_path(color="white") +
  coord_equal() +
  scale_fill_brewer(name = NULL,
                    type="qual",
                    palette = "Set3") +
  geom_polygon(data = garden_locations.df, mapping = aes(long,lat,group=group), inherit.aes = F) +
  geom_text(data = garden_locations %>%
              rgeos::gCentroid(byid=T) %>% as.data.frame() %>%
              cbind(garden_locations@data),
            mapping = aes(x,
                          y,
                          label = Abbreviation,
                          vjust = 0,
                          hjust = c(0,0,1,0)), 
            inherit.aes = F,
            nudge_x = c(10,10,-10,10),
            nudge_y = 10,
            size=3, color = "black"
              ) +
  xlab("Easting") +
  ylab("Northing") +
  xlim(710000,711600) +
  ylim(4136400,4137800) +
  theme(axis.text=element_text(size=8, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))

##### END garden_soils_map.pdf #####

##### figure_4.pdf #####
## A grid of maps showing the salient soil characteristics
## Maps of soil characteristics included in the Village Ecodynamics Project 
## maize paleoproductivity reconstructions.
# Prepare the data for plotting in ggplot2
ccac_nad27 <- ccac %>% sp::spTransform(sp::CRS(raster::projection(vepi.soils.ccac.raster)))

vep_soil_raster_soils <- rasterVis::gplot(vepi.soils.ccac.raster[["VPSCode"]] %>% raster::crop(ccac_nad27)) + 
  geom_raster(aes(fill=factor(value, levels = sort(unique(soils_VEPI_yields.df$VPSCode)))), na.rm = T) +
  ggtitle("Rasterized soils") +
  coord_equal() +
  scale_fill_brewer(name = NULL,
                    type="qual",
                    palette = "Set3",
                    na.translate = FALSE,
                    drop = FALSE) +
  geom_polygon(data = garden_locations.df, mapping = aes(long,lat,group=group), inherit.aes = F) +
  geom_text(data = garden_locations %>%
              rgeos::gCentroid(byid=T) %>% as.data.frame() %>%
              cbind(garden_locations@data),
            mapping = aes(x,
                          y,
                          label = Abbreviation,
                          vjust = 0,
                          hjust = c(0,0,1,0)), 
            inherit.aes = F,
            nudge_x = c(10,10,-10,10),
            nudge_y = 10,
            size=3, color = "black"
  ) +
  xlab("Easting") +
  ylab("Northing") +
  xlim(710000,711600) +
  ylim(4136400,4137800) +
  theme(axis.text=element_text(size=8, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))

vep_soil_raster_clusters <- rasterVis::gplot(vepi.soils.ccac.raster[["CLUSTER"]]) + 
  geom_raster(aes(fill=factor(value, levels = sort(unique(value)))), na.rm = T) +
  ggtitle("VEP soil clusters") +
  coord_equal() +
  scale_fill_brewer(name = NULL,
                    type="qual",
                    palette = "Set3",
                    na.translate = FALSE,
                    drop = FALSE) +
  geom_polygon(data = garden_locations.df, mapping = aes(long,lat,group=group), inherit.aes = F) +
  geom_text(data = garden_locations %>%
              rgeos::gCentroid(byid=T) %>% as.data.frame() %>%
              cbind(garden_locations@data),
            mapping = aes(x,
                          y,
                          label = Abbreviation,
                          vjust = 0,
                          hjust = c(0,0,1,0)), 
            inherit.aes = F,
            nudge_x = c(10,10,-10,10),
            nudge_y = 10,
            size=3, color = "black"
  ) +
  xlab("Easting") +
  ylab("Northing") +
  xlim(710000,711600) +
  ylim(4136400,4137800) +
  theme(axis.text=element_text(size=8, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))

multiplot <- function(vars, limits, titles){
  mapply(vars, limits, titles, SIMPLIFY = F, FUN = function(the.var,the.limit,the.title){
    out <- rasterVis::gplot(vepi.soils.ccac.raster[[the.var]]) + 
      geom_raster(aes(fill = value), na.rm = T) +
      ggtitle(the.title) +
      coord_equal() +
      scale_fill_distiller(name = NULL,
                           limits = the.limit,
                           palette = "YlGn",
                           direction = 1,
                           na.value = NA) +
      geom_polygon(data = garden_locations.df, mapping = aes(long,lat,group=group), inherit.aes = F) +
      geom_text(data = garden_locations %>%
                  rgeos::gCentroid(byid=T) %>% as.data.frame() %>%
                  cbind(garden_locations@data),
                mapping = aes(x,
                              y,
                              label = Abbreviation,
                              vjust = 0,
                              hjust = c(0,0,1,0)), 
                inherit.aes = F,
                nudge_x = c(10,10,-10,10),
                nudge_y = 10,
                size=3, color = "black"
      ) +
      xlab("Easting") +
      ylab("Northing") +
      xlim(710000,711600) +
      ylim(4136400,4137800) +
      theme(axis.text=element_text(size=8, color = "black"),
            title=element_text(size=8,face="bold", color = "black"),
            legend.text=element_text(size=6, color = "black"))
    return(ggplotGrob(out))
  })
}
vars <- c("NYProd_lb_ac","AWC_Lower_median","SCM_RED")
limits <- list(c(0,1250), c(0,10), c(0,1))
titles <- c("Normal-year productivity (lb/ac)", "AWC 6–60 inches (in)", "Hand planting factor")

the.plots <- multiplot(vars, limits, titles)

the.plots <- c(list(vep_soil_polygons = ggplotGrob(vep_soil_polygons),vep_soil_raster_soils = ggplotGrob(vep_soil_raster_soils),vep_soil_raster_clusters = ggplotGrob(vep_soil_raster_clusters)),the.plots)

mai <- c(0,0,0,0)
fig.width <- 7.5
fig.height <- 8.5

cairo_pdf(file="./OUTPUT/FIGURES/figure_4.pdf",
          width=fig.width,
          height=fig.height,
          antialias="none",
          bg="white",
          pointsize=6)
par(mai=mai, xpd=F)
gridExtra::grid.arrange(grobs = the.plots, nrow=3)
dev.off()
distill("./OUTPUT/FIGURES/figure_4.pdf")
##### END figure_4.pdf #####


##### table_2.csv #####
## Basic data on each garden
garden_yields <- yields %>%
  dplyr::select(Garden, Season, 
                # `Yield by clump area`, 
                `PFP experimental yield (kg/ha)`) %>%
  dplyr::group_by(Garden, Season) %>%
  dplyr::summarise(
    # `Yield by clump area` = mean(`Yield by clump area`), 
    `PFP experimental yield (kg/ha)` = mean(`PFP experimental yield (kg/ha)`)
  ) %>%
  dplyr::left_join(yields %>%
                     dplyr::select(Garden, Season, Variety, Clumps, `Spacing (m)`) %>%
                     unique(), by = c("Garden","Season")) %>%
  dplyr::select(Garden, Season, Variety, Clumps, `Spacing (m)`,
                # `Yield by clump area`,
                `PFP experimental yield (kg/ha)`) %>%
  dplyr::left_join(PFP_VEP_yields, by=c("Garden","Season")) %>%
  dplyr::rename(`VEP estimated yield (kg/ha)` = Yield)

garden_yields %>% 
    dplyr::mutate(`PFP experimental yield (kg/ha)` = round(`PFP experimental yield (kg/ha)`, digits = 1),
                  `VEP estimated yield (kg/ha)` = round(`VEP estimated yield (kg/ha)`, digits = 1)) %>%
      write_csv(path = "./OUTPUT/TABLES/table_2.csv")

##### END table_2.csv #####

##### table_1.csv #####
soils_data <- garden_locations %>%
  rgeos::gCentroid(byid=T) %>% # Get the centroid under each garden
  sp::spTransform(sp::CRS(raster::projection(vepi.soils.ccac.raster))) %>% # transform into the same projection
  raster::extract(x = vepi.soils.ccac.raster) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(Garden = garden_locations$Abbreviation) %>%
  dplyr::left_join(vepi.soils.ccac@data %>%
                     dplyr::select(VPSCode, muname, MUKEY) %>%
                     dplyr::distinct(), by = "VPSCode")

soils_data %>%
  dplyr::select(muname, MUKEY, VPSCode, CLUSTER, Garden, AWC_Lower_median:SCM_RED) %>%
  dplyr::rename(`Map unit name` = muname,
                `Map unit key` = MUKEY,
                `VEP soil ID` = VPSCode,
                `VEP soil cluster` = CLUSTER,
                `Normal-year prod. (lb/ac)` = NYProd_lb_ac,
                `AWC 6–60 inches (in)` = AWC_Lower_median,
                `Hand planting factor` = SCM_RED,
                `PFP Garden` = Garden) %>%
  dplyr::mutate(`AWC 6–60 inches (in)` = round(`AWC 6–60 inches (in)`, digits = 2),
                `Hand planting factor` = round(`Hand planting factor`, digits = 2),
                `Normal-year prod. (lb/ac)` = round(`Normal-year prod. (lb/ac)`, digits = 1)) %>%
  write_csv(path = "./OUTPUT/TABLES/table_1.csv")
##### END table_1.csv #####

##### figure_5.pdf #####
## A visual comparison between the experimental and VEP estimated yields
## Experimental (PFP) and estimated (VEP) garden yields.
## Box plots indicate the distribution of experimental yields as extrapolated from individual clumps.
## Asterisks mark the distribution means.
## The CDG, POG, and PLC gardens all occur in the same NRCS soil mapunit, 
## and therefor share the same estimated yield.
mai <- c(0.25,0.25,0,0)
fig.width <- 5
fig.height <- 2.25

cairo_pdf(file="./OUTPUT/FIGURES/figure_5.pdf", width=fig.width, height=fig.height, antialias="none", bg="white", pointsize=6, fallback_resolution = 600)
par(mai=mai, xpd=F)

ggplot() + 
  geom_boxplot(data = yields,
               mapping = aes(y = `PFP experimental yield (kg/ha)`, x = Season, fill = Garden),
               outlier.size = 0.5,
               size = 0.25,
               colour = "black") +
  scale_fill_brewer("Garden Yields:\nPFP Experiments", palette="OrRd") +
  stat_summary(data = yields,
               mapping = aes(y = `PFP experimental yield (kg/ha)`, x = Season, fill = Garden),
               fun.y=mean,
               colour="black", 
               geom="point",
               position=position_dodge(width=0.75),
               size = 0.75,
               shape = 8) +
  geom_line(data = PFP_VEP_yields,
            mapping = aes(x=Season, y=Yield, group = Garden, colour = Garden), 
            size = 0.25
            ) +
  scale_colour_brewer("Garden Yields:\nVEP Estimates", palette="OrRd") +
  ylab("Yield (kg/ha)") +
  theme(axis.text=element_text(size=6, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))

dev.off()
distill("./OUTPUT/FIGURES/figure_5.pdf")

##### END figure_5.pdf #####

## Calculate correlation between each garden's 
## mean experimental and estimated yield
yield_correlations <- garden_yields %>% 
  do(tidy(cor.test(.$`VEP estimated yield (kg/ha)`, .$`PFP experimental yield (kg/ha)`))) %>%
  dplyr::select(Garden, estimate, p.value, conf.low, conf.high) %>%
  dplyr::rename(`Correlation (r)` = estimate, 
                `P-value` = p.value, 
                `Lower CI` = conf.low, 
                `Upper CI` = conf.high) 
yield_correlations %>% 
  dplyr::mutate(`Correlation (r)` = round(`Correlation (r)`, digits = 2),
                `P-value` = round(`P-value`, digits = 3),
                `Lower CI` = round(`Lower CI`, digits = 2),
                `Upper CI` = round(`Upper CI`, digits = 2)) %>%
  write_csv(path = "./OUTPUT/TABLES/table_4.csv")

##### Yield mean variance table #####
garden_prod_vepi <- garden_locations %>%
  rgeos::gCentroid(byid=T) %>% # Get the centroid under each garden
  sp::spTransform(sp::CRS(raster::projection(vepi.paleoprod))) %>% # transform into the same projection
  raster::extract(x = vepi.paleoprod) %>%
  t() %>%
  tibble::as_tibble()  %>%
  set_colnames(garden_locations$Abbreviation) %>%
  dplyr::mutate(Year = 600:1300) %>% 
  tidyr::gather(Garden,Yield,CDG:PLC) %>%
  dplyr::group_by(Garden) %>%
  dplyr::summarise(`VEP estimated yield: mean, ancient` = mean(`Yield`),
                   `VEP estimated yield: SD, ancient` = sd(`Yield`)
  )

yield_mean_variance <- garden_yields %>%
  dplyr::summarise_each(funs = c("mean","sd"), `PFP experimental yield (kg/ha)`,`VEP estimated yield (kg/ha)`) %>%
  dplyr::rename(`PFP experimental yield: mean` = `PFP experimental yield (kg/ha)_mean`, 
                `VEP estimated yield: mean, modern` = `VEP estimated yield (kg/ha)_mean`,
                `PFP experimental yield: SD` = `PFP experimental yield (kg/ha)_sd`, 
                `VEP estimated yield: SD, modern` = `VEP estimated yield (kg/ha)_sd`) %>%
  dplyr::left_join(garden_prod_vepi, by = "Garden") %>%
  dplyr::select(Garden,
                `PFP experimental yield: mean`,
                `VEP estimated yield: mean, modern`,
                `VEP estimated yield: mean, ancient`,
                `PFP experimental yield: SD`,
                `VEP estimated yield: SD, modern`,
                `VEP estimated yield: SD, ancient`) 
yield_mean_variance %>%
  dplyr::mutate(Garden,
                `PFP experimental yield: mean` = round(`PFP experimental yield: mean`, digits = 1),
                `VEP estimated yield: mean, modern` = round(`VEP estimated yield: mean, modern`, digits = 1),
                `VEP estimated yield: mean, ancient` = round(`VEP estimated yield: mean, ancient`, digits = 1),
                `PFP experimental yield: SD` = round(`PFP experimental yield: SD`, digits = 1),
                `VEP estimated yield: SD, modern` = round(`VEP estimated yield: SD, modern`, digits = 1),
                `VEP estimated yield: SD, ancient` = round(`VEP estimated yield: SD, ancient`, digits = 1)) %>%
  dplyr::rename(`PFP experimental yield: mean (kg/ha)` = `PFP experimental yield: mean`,
                `VEP estimated yield: mean, modern (kg/ha)` = `VEP estimated yield: mean, modern`,
                `VEP estimated yield: mean, ancient (kg/ha)` = `VEP estimated yield: mean, ancient`,
                `PFP experimental yield: SD (kg/ha)` = `PFP experimental yield: SD`,
                `VEP estimated yield: SD, modern (kg/ha)` = `VEP estimated yield: SD, modern`,
                `VEP estimated yield: SD, ancient (kg/ha)` = `VEP estimated yield: SD, ancient`) %>%
  write_csv(path = "./OUTPUT/TABLES/table_3.csv")

##### figure_3.pdf #####
# Force Raster to load large rasters into memory
raster::rasterOptions(chunksize=2e+08,maxmemory=2e+09)
vepi.paleoprod <- vepi.paleoprod * 1

space_time_plot(the_brick = vepi.paleoprod,
                out_file = "./OUTPUT/FIGURES/figure_3.pdf",
                timelim = c(600,1300),
                timeaxis = seq(700,1300,100),
                zlim_mid_range = c(0,500),
                zlab = "Yeild (kg/ha)",
                zaxis = c(100,200,300,400))
##### END figure_3.pdf #####

# Distill Burns and Van West figures.
file.copy("./DATA/BURNS.pdf","./OUTPUT/FIGURES/figure_1.pdf", overwrite = TRUE)
distill("./OUTPUT/FIGURES/figure_1.pdf")
file.copy("./DATA/VANWEST.pdf","./OUTPUT/FIGURES/figure_2.pdf", overwrite = TRUE)
distill("./OUTPUT/FIGURES/figure_2.pdf")

##### Zip up SI data #####
# zip("./Bocinsky_Varien_2017.zip", c("./DATA/","./Bocinsky_Varien_2017.R","./Bocinsky_Varien_2017.Rproj","./src/","./README.md"), flags = "-r9X", extras = "-r --exclude=*.DS_Store* --exclude=*.git*")
##### END Zip up SI data #####
