## This is the script for the analyses in Varien and Bocinsky, under review:
## Varien, Mark D. and R. Kyle Bocinsky. Under review. Calibrating Maize 
## Paleoproduction Models using Experimental Data. Ethnobiology.

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
## Date: 8/31/2016

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
FedData::pkg_test("ggplot2")
FedData::pkg_test("gridExtra")
FedData::pkg_test("raster")
FedData::pkg_test("sp")
FedData::pkg_test("RColorBrewer")
FedData::pkg_test("rgeos")
FedData::pkg_test("broom")
FedData::pkg_test("maptools")
FedData::pkg_test("dplyr")

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
  unlink("./Varien_Bocinsky_2017.zip", force = T)
}

# Create directories for output
dir.create("./OUTPUT/DATA", recursive = T, showWarnings = F)
dir.create("./OUTPUT/FIGURES", recursive = T, showWarnings = F)
dir.create("./OUTPUT/TABLES", recursive = T, showWarnings = F)

# Growing seasons of the PFP
seasons <- 2009:2015


##### BEGIN VEP PRODUCTIVITY ESTIMATION #####

### SOILS ###

# Load a polygon of the Indian Camp Ranch + Crow Canyon study area
CCAC <- rgdal::readOGR(dsn = "./DATA/ccac.geojson", "OGRGeoJSON", verbose = FALSE)

# Download the NRCS SSURGO soils data for the study area
# This downloads the soils data for the CO671 soil survey
CCAC_SSURGO <- FedData::get_ssurgo(template=CCAC, label="CCAC", raw.dir="./OUTPUT/DATA/SSURGO/RAW/", extraction.dir="./OUTPUT/DATA/SSURGO/EXTRACTIONS/", force.redo=F)

## Convert all tables to tibbles
CCAC_SSURGO$tabular <- lapply(CCAC_SSURGO$tabular,dplyr::as_data_frame)

# Calculate all mapunit-level NPP, bean productivity, available water content, and hand planting restriction
CCAC_SSURGO <- CCAC_SSURGO %>%
  calculate_mukey_annual_npp() %>% # Calculate annual NPP
  calculate_mukey_bean_productivity() %>% # Calculate bean productivity
  calculate_mukey_awc() %>% # Calculate and output AWC
  calculate_mukey_hand_planting_restriction() # Calculate hand-planting restriction

# Add the data to the shapefile of soil mapunits
soils <- CCAC_SSURGO$spatial
soils@data <- dplyr::left_join(soils@data,CCAC_SSURGO$tabular$mapunit %>%
                                 dplyr::select(mukey,muname) %>%
                                 dplyr::mutate(mukey = as.factor(mukey)),
                               by = c("MUKEY" = "mukey"))

# Calculate area of each mapunit
soils$area <- rgeos::gArea(spTransform(soils,CRS("+proj=utm +datum=NAD83 +zone=12")),byid=T)

# Clean up the data
soils@data <- soils@data %>%
  dplyr::select(-area,-AREASYMBOL,-SPATIALVER, -MUSYM) %>%
  dplyr::as_data_frame() %>%
  as.data.frame()

unlink('./OUTPUT/DATA/soils')
rgdal::writeOGR(soils,dsn='./OUTPUT/DATA/soils',layer="soils",driver="GeoJSON",overwrite=T)
file.rename('./OUTPUT/DATA/soils','./OUTPUT/DATA/soils.geojson')

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
                                                    raw.dir = "./OUTPUT/GHCN"),
                    FedData::get_ghcn_daily_station(ID="USC00051886", 
                                                    elements = c("PRCP"), 
                                                    raw.dir = "./OUTPUT/GHCN"))

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

## Read in the soils data generated above
soils <- rgdal::readOGR(dsn = "./OUTPUT/DATA/soils.geojson", "OGRGeoJSON", verbose = FALSE)

# Calculate the centroid of each soil and add to the data table
soils_data <- cbind(soils@data, as.data.frame(rgeos::gCentroid(soils, byid=T)))

## Produce PDSI estimates for each soil
# Create an output directory for PDSI data
dir.create("./OUTPUT/DATA/PDSI/", showWarnings = FALSE, recursive = TRUE)
output.dir <- "./OUTPUT/DATA/PDSI/"
scPDSI.path <- "./src/scpdsi"

# Create a "blank" PDSI brick
# Calculate June Monthly PDSI for each soil
# THIS WILL ONLY RUN ON A MAC OR OTHER UNIX-ALIKE
monthly_T <- data.frame(year=unique(cortez_weather_monthly$YEAR),matrix(paste0(" ",gsub("0000NA","",sprintf("%06.3f",cortez_weather_monthly$TAVG_F))),ncol=12,byrow=T), stringsAsFactors=F)
monthly_P <- data.frame(year=unique(cortez_weather_monthly$YEAR),matrix(paste0(" ",gsub("0000NA","",sprintf("%06.3f",cortez_weather_monthly$PRCP_IN))),ncol=12,byrow=T), stringsAsFactors=F)
mon_T_normal <- data.frame(matrix(paste0(" ",as.character(formatC(cortez_weather_monthly_norms$TAVG_F, format = 'f', digits = 3, width=6))),ncol=12,byrow=T), stringsAsFactors=F)
mon_T_normal[1,1] <- paste0(" ",mon_T_normal[1,1])
Monthly_PDSI <- lapply(1:nrow(soils_data),FUN = function(x){
  return(rPDSI(output.dir = output.dir, monthly_T = monthly_T, monthly_P = monthly_P, mon_T_normal = mon_T_normal, awc = soils_data[x,'lower_awc'], lat = soils_data[x,'y'], scPDSI.path = scPDSI.path))
}) %>%
  lapply(FUN = function(x){
    out <- dplyr::bind_cols(cortez_weather_monthly %>% dplyr::select(YEAR,MONTH),x) %>%
      dplyr::filter(YEAR %in% seasons, MONTH == 6) %>%
      dplyr::select(YEAR,PDSI) %>%
      collect %>%
      .[["PDSI"]]
    names(out) <- seasons
    return(out)
  })


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
predictions <- lapply(Monthly_PDSI,function(x){
  (pdsi.coef*x + intercept) * 62.77 # Convert to kg/ha
})
predictions <- do.call(rbind,predictions)

## Re-weight production for each soil 
# Read in the NRCS normal year dry-weight soil productivity (NPP) & Bean Soils
NPP.bean.mean <- 1093 * 1.12085 # Mean Bean Soil NPP for all VEP I soils, converted to kg/ha 
NPP.reweight <- soils_data$NPP/NPP.bean.mean
predictions <- sweep(predictions,MARGIN=1,NPP.reweight,`*`)

## We depart a bit from the methods described in Kohler 2012. Before, the retrodicted values
## would be multiplied by the NPP renormed values, then loaded into the simulation, where they 
## would be renormed to prehispanic maize varieties and further reduced by a cold correction 
## and a hand-planting reduction (steps 9, 10, and 11 in Kohler 2012).
## We perform all of those steps "up front" here.

## RENORM MAIZE PRODUCTION FOR PREHISPANIC VARIETIES AND CULTIVATION PRACTICES
# Kohler 2012: Step 9 (pp. 100--103)
renorm.factor <- 0.68
predictions <- predictions * 0.68

## MAKE HAND-PLANTING ADJUSTMENT
# Kohler 2012: Step 10 (pp. 103)
# Read in hand-planting factor from the soils data
predictions <- sweep(predictions,MARGIN=1,soils_data$SCM_RED,`*`)


## MAKE COLD CORRECTION
# Kohler 2012: Step 11 (pp. 103--106)
# We don't make the cold correction here, because CCAC PFP gardens are below 7,000 feet.

# Write yield info to soils shapefile
soils@data <- cbind(soils@data,predictions)

unlink('./OUTPUT/DATA/soils_VEPII_yields')
rgdal::writeOGR(soils,dsn='./OUTPUT/DATA/soils_VEPII_yields',layer="soils_VEPII_yields",driver="GeoJSON",overwrite=T)
file.rename('./OUTPUT/DATA/soils_VEPII_yields','./OUTPUT/DATA/soils_VEPII_yields.geojson')

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
  dplyr::mutate(`Yield by clump area` = (Clumps * `Net kernel weight`/1000)/((Spacing ^ 2) * Clumps * 0.0001), # yield by actual clump area
                `Standardized yield` = (Clumps * `Net kernel weight`/1000)/((2 ^ 2) * Clumps * 0.0001) # yield by clump area with 2m spacing
                ) %>% 
  dplyr::mutate(Variety = as.factor(Variety),
                Garden = as.factor(Garden),
                Season = as.factor(Season)) %>%
  dplyr::arrange(Season, Garden, Clump) %>% # sort by these variables
  dplyr::select(Season, Garden, Variety, Clumps, Spacing, Clump, `Net kernel weight`, `Yield by clump area`, `Standardized yield`) %>% # reorder columns
  dplyr::group_by(Season, Garden) # calculations are by season and garden

readr::write_csv(yields,"./OUTPUT/DATA/yields.csv")

##### END PFP PRODUCTIVITY ESTIMATION #####

##### BEGIN JOINT ANALYSIS #####
## Extract the VEP yield reconstructions under the PFP gardens
# Read in the VEP productivity reconstruction from above
soils_VEPII_yields <- rgdal::readOGR(dsn = "./OUTPUT/DATA/soils_VEPII_yields.geojson", "OGRGeoJSON", verbose = FALSE)
soils_VEPII_yields$muname <- gsub(" MLRA 36","",soils_VEPII_yields$muname)

PFP_VEP_yields <- garden_locations %>%
  rgeos::gCentroid(byid=T) %>% # Get the centroid under each garden
  sp::spTransform(sp::CRS(raster::projection(soils_VEPII_yields))) %>% # transform into the same projection
  raster::extract(x = soils_VEPII_yields) %>%
  dplyr::mutate(Garden = garden_locations$Abbreviation) %>%
  dplyr::select(Garden,MUKEY:SCM_RED,X2009:X2015) %>%
  tidyr::gather(Season, Yield, X2009:X2015) %>%
  dplyr::as_data_frame() %>%
  dplyr::mutate(Season = gsub("X","",Season) %>% as.factor()) %>%
  dplyr::select(Garden,Season,Yield) %>%
  dplyr::arrange(Garden, Season)

##### garden_soils_map.pdf #####
## A map showing the soils on the CCAC campus and the locations of the PFP gardens
# Prepare the data for plotting in ggplot2
soils_VEPII_yields@data$id = rownames(soils_VEPII_yields@data)
soils_VEPII_yields.points = ggplot2::fortify(soils_VEPII_yields %>% sp::spTransform(sp::CRS("+proj=utm +nad=NAD83 +zone=12")), region="id") %>% as_data_frame()
soils_VEPII_yields.df = dplyr::full_join(soils_VEPII_yields.points, soils_VEPII_yields@data, by="id")
garden_locations@data$id = rownames(garden_locations@data)
garden_locations.points = ggplot2::fortify(garden_locations %>% sp::spTransform(sp::CRS("+proj=utm +nad=NAD83 +zone=12")), region="id") %>% as_data_frame()
garden_locations.df = dplyr::full_join(garden_locations.points, garden_locations@data, by="id")

mai <- c(0,0,0,0)
fig.width <- 6.5
fig.height <- 4

quartz(file="./OUTPUT/FIGURES/garden_soils_map.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=6, dpi=600)
par(mai=mai, xpd=F)
ggplot(soils_VEPII_yields.df) +
  aes(long,lat,group=group,fill=paste0(gsub(", ",",\n",muname), "\n")) +
  geom_polygon() + 
  geom_path(color="white") +
  coord_equal() +
  scale_fill_brewer(name = NULL, type="qual", palette = "Set3") +
  geom_polygon(data = garden_locations.df, mapping = aes(long,lat,group=group), inherit.aes = F) +
  geom_text(data = garden_locations %>%
              sp::spTransform(sp::CRS("+proj=utm +nad=NAD83 +zone=12")) %>%
              rgeos::gCentroid(byid=T) %>% as.data.frame() %>%
              cbind(garden_locations@data),
            mapping = aes(x,
                          y,
                          label = Abbreviation,
                          vjust = 0,
                          hjust = c(0,0,1,0)), 
            inherit.aes = F,
            nudge_x = c(10,10,-10,10),
            nudge_y = 10
              ) +
  xlab("Easting") +
  ylab("Northing") +
  theme(axis.text=element_text(size=8, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))
dev.off()
distill("./OUTPUT/FIGURES/garden_soils_map.pdf")

##### END garden_soils_map.pdf #####

##### soils_vars_grid.pdf #####
## A grid of maps showing the salient soil characteristics
# Prepare the data for plotting in ggplot2
soils_VEPII_yields@data$id = rownames(soils_VEPII_yields@data)
soils_VEPII_yields.points = ggplot2::fortify(soils_VEPII_yields %>% sp::spTransform(sp::CRS("+proj=utm +nad=NAD83 +zone=12")), region="id") %>% as_data_frame()
soils_VEPII_yields.df = dplyr::full_join(soils_VEPII_yields.points, soils_VEPII_yields@data, by="id")

multiplot <- function(vars, limits, titles){
  mapply(vars, limits, titles, SIMPLIFY = F, FUN = function(the.var,the.limit,the.title){
    out <- ggplot(soils_VEPII_yields.df) + 
      ggtitle(the.title) +
      aes_string("long","lat",group="group",fill=the.var) + 
      geom_polygon() +
      geom_path(color="white") +
      coord_equal() +
      scale_fill_distiller(name = NULL, limits = the.limit, palette = "YlGn", direction = 1) +
      xlab("Easting") +
      ylab("Northing") +
      theme(axis.text=element_text(size=8, color = "black"),
            title=element_text(size=8,face="bold", color = "black"),
            legend.text=element_text(size=6, color = "black"))
    return(ggplotGrob(out))
  })
}
vars <- c("NPP","bean_yield","lower_awc","SCM_RED")
limits <- list(c(0,1250), c(0,500), c(0,8), c(0.85,1))
titles <- c("Net primary productivity (lb/ac)", "Bean yield (lb/ac)", "AWC 6–60 inches (in/in)", "Hand planting factor")

the.plots <- multiplot(vars, limits, titles)

mai <- c(0,0,0,0)
fig.width <- 9
fig.height <- 6.5

quartz(file="./OUTPUT/FIGURES/soils_vars_grid.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=6, dpi=600)
par(mai=mai, xpd=F)
gridExtra::grid.arrange(grobs = the.plots, nrow=2)
dev.off()
distill("./OUTPUT/FIGURES/soils_vars_grid.pdf")

##### END soils_vars_grid.pdf #####


##### garden_yields.csv #####
## Basic data on each garden
yields %>%
  dplyr::select(Garden, Season, `Yield by clump area`, `Standardized yield`) %>%
  dplyr::group_by(Garden, Season) %>%
  dplyr::summarise(`Yield by clump area` = mean(`Yield by clump area`), `Standardized yield` = mean(`Standardized yield`)) %>%
  dplyr::left_join(yields %>%
                     dplyr::select(Garden, Season, Variety, Clumps, Spacing) %>%
                     unique(), by = c("Garden","Season")) %>%
  dplyr::select(Garden, Season, Variety, Clumps, Spacing, `Yield by clump area`, `Standardized yield`) %>%
  dplyr::left_join(PFP_VEP_yields, by=c("Garden","Season")) %>%
  dplyr::rename(`VEP yield estimate` = Yield) %>%
  write_csv(path = "./OUTPUT/TABLES/garden_yields.csv")
##### END garden_yields.csv #####

##### clump_yields.csv #####
## Basic data on each clump
yields %>%
  dplyr::select(Garden, Season, Clump, `Net kernel weight`) %>%
  dplyr::arrange(Garden, Season, Clump) %>%
  write_csv(path = "./OUTPUT/TABLES/clump_yields.csv")
##### END clump_yields.csv #####

##### soils_data.csv #####
soils_data <- garden_locations %>%
  rgeos::gCentroid(byid=T) %>% # Get the centroid under each garden
  sp::spTransform(sp::CRS(raster::projection(soils_VEPII_yields))) %>% # transform into the same projection
  raster::extract(x = soils_VEPII_yields) %>%
  dplyr::mutate(Garden = garden_locations$Abbreviation)

soils_data %>%
  dplyr::select(MUKEY, Garden) %>%
  dplyr::group_by(MUKEY) %>%
  dplyr::summarise(Garden = paste0(Garden, collapse = ", ")) %>%
  dplyr::left_join(soils_data %>% dplyr::select(-Garden, -point.ID, -poly.ID, -upper_awc), by = ("MUKEY")) %>%
  unique() %>%
  dplyr::select(muname, MUKEY, Garden, NPP:SCM_RED) %>%
  dplyr::rename(`Map unit name` = muname, `Map unit key` = MUKEY, `Net primary prod. (lb/ac)` = NPP, `Bean yield (lb/ac)` = bean_yield, `AWC 6–60 inches (in/in)` = lower_awc, `Hand planting factor` = SCM_RED, `PFP Gardens` = Garden) %>%
  write_csv(path = "./OUTPUT/TABLES/soils_data.csv")
##### END soils_data.csv #####

##### yields_standardized.pdf #####
## A visual comparison between the experimental and VEP estimated yields
mai <- c(0.25,0.25,0,0)
fig.width <- 5
fig.height <- 2

quartz(file="./OUTPUT/FIGURES/yields_standardized.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=6, dpi=600)
par(mai=mai, xpd=F)

ggplot() + 
  geom_boxplot(data = yields, mapping = aes(y = `Standardized yield`, x = Season, fill = Garden), outlier.size = 0.5, size = 0.25, colour = "black") +
  scale_fill_brewer("Garden Yields:\nExperimental", palette="OrRd") +
  stat_summary(data = yields,
               mapping = aes(y = `Standardized yield`, x = Season, fill = Garden),
               fun.y=mean,
               colour="black", 
               geom="point",
               position=position_dodge(width=0.75),
               size = 0.75,
               shape = 8) +
  geom_line(data = PFP_VEP_yields %>%
              dplyr::mutate(Garden = ifelse(Garden == "KUG","KUG","CDG, POG, PLC")) %>%
              unique(),
            mapping = aes(x=Season, y=Yield, group = Garden, linetype = Garden), 
            size = 0.25
            ) +
  scale_linetype("Garden Yields:\nVEP Estimates") +
  ylab("Yield (kg/ha)") +
  theme(axis.text=element_text(size=6, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))

dev.off()
distill("./OUTPUT/FIGURES/yields_standardized.pdf")

##### END yields_standardized.pdf #####

##### yields_clump_spacing.pdf #####
## A visual comparison between the experimental and VEP estimated yields
mai <- c(0.25,0.25,0,0)
fig.width <- 5
fig.height <- 2

quartz(file="./OUTPUT/FIGURES/yields_clump_spacing.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=6, dpi=600)
par(mai=mai, xpd=F)

ggplot() + 
  geom_boxplot(data = yields, mapping = aes(y = `Yield by clump area`, x = Season, fill = Garden), outlier.size = 0.5, size = 0.25, colour = "black") +
  scale_fill_brewer("Garden Yields:\nExperimental", palette="OrRd") +
  stat_summary(data = yields,
               mapping = aes(y = `Yield by clump area`, x = Season, fill = Garden),
               fun.y=mean,
               colour="black", 
               geom="point",
               position=position_dodge(width=0.75),
               size = 0.75,
               shape = 8) +
  geom_line(data = PFP_VEP_yields %>%
              dplyr::mutate(Garden = ifelse(Garden == "KUG","KUG","CDG, POG, PLC")) %>%
              unique(),
            mapping = aes(x=Season, y=Yield, group = Garden, linetype = Garden), 
            size = 0.25
  ) +
  scale_linetype("Garden Yields:\nVEP Estimates") +
  ylab("Yield (kg/ha)") +
  theme(axis.text=element_text(size=6, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))

dev.off()
distill("./OUTPUT/FIGURES/yields_clump_spacing.pdf")

##### END yields_clump_spacing.pdf #####

## Calculate correlation between each garden's 
## mean experimental and estimated yield
yields %>%
  dplyr::select(Garden, Season, `Standardized yield`) %>%
  dplyr::group_by(Garden, Season) %>%
  dplyr::summarise(`Standardized yield` = mean(`Standardized yield`)) %>%
  dplyr::full_join(PFP_VEP_yields, by = c("Garden","Season")) %>%
  dplyr::rename(`Experimental yield` = `Standardized yield`, `Estimated yield` = Yield) %>%
  dplyr::mutate(`Experimental yield` = scale(`Experimental yield`) %>% as.numeric(), 
                `Estimated yield` = scale(`Estimated yield`) %>% as.numeric()) %>%
  dplyr::group_by(Garden) %>% 
  do(tidy(cor.test(.$`Estimated yield`, .$`Experimental yield`))) %>%
  dplyr::select(Garden, estimate, p.value, conf.low, conf.high) %>%
  dplyr::rename(Correlation = estimate, 
                `P-value` = p.value, 
                `Lower CI` = conf.low, 
                `Upper CI` = conf.high) %T>%
  write_csv(path = "./OUTPUT/TABLES/yield_correlations_standardized.csv") %>%
  dplyr::ungroup() %>%
  dplyr::summarise(mean(Correlation))

##### yield_correlations_standardized.pdf #####
## A visual comparison between the scaled experimental and VEP estimated yields
## with linear regression
phi <- (1+sqrt(5))/2
mai <- c(0.25,0.25,0,0)
fig.width <- 5
fig.height <- fig.width/phi

quartz(file="./OUTPUT/FIGURES/yield_correlations_standardized.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=6, dpi=600)
par(mai=mai, xpd=F)

yields %>%
  dplyr::select(Garden, Season, `Standardized yield`) %>%
  dplyr::group_by(Garden, Season) %>%
  dplyr::summarise(`Standardized yield` = mean(`Standardized yield`)) %>%
  dplyr::full_join(PFP_VEP_yields, by = c("Garden","Season")) %>%
  dplyr::rename(`Experimental yield` = `Standardized yield`, `Estimated yield` = Yield) %>%
  dplyr::mutate(`Experimental yield` = scale(`Experimental yield`) %>% as.numeric(), 
                `Estimated yield` = scale(`Estimated yield`) %>% as.numeric()) %>%
  ggplot(aes(x = `Estimated yield`, y = `Experimental yield`, color = Garden)) + 
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE,
              size = 0.5) +   # Don't add shaded confidence region
  geom_point() +
  scale_colour_brewer("Garden", palette="OrRd") +
  theme(axis.text=element_text(size=6, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))

dev.off()
distill("./OUTPUT/FIGURES/yield_correlations_standardized.pdf")

##### END yield_correlations_standardized.pdf #####

## Calculate correlation between each garden's 
## mean experimental and estimated yield
yields %>%
  dplyr::select(Garden, Season, `Yield by clump area`) %>%
  dplyr::group_by(Garden, Season) %>%
  dplyr::summarise(`Yield by clump area` = mean(`Yield by clump area`)) %>%
  dplyr::full_join(PFP_VEP_yields, by = c("Garden","Season")) %>%
  dplyr::rename(`Experimental yield` = `Yield by clump area`, `Estimated yield` = Yield) %>%
  dplyr::mutate(`Experimental yield` = scale(`Experimental yield`) %>% as.numeric(), 
                `Estimated yield` = scale(`Estimated yield`) %>% as.numeric()) %>%
  dplyr::group_by(Garden) %>% 
  do(tidy(cor.test(.$`Estimated yield`, .$`Experimental yield`))) %>%
  dplyr::select(Garden, estimate, p.value, conf.low, conf.high) %>%
  dplyr::rename(Correlation = estimate, 
                `P-value` = p.value, 
                `Lower CI` = conf.low, 
                `Upper CI` = conf.high) %T>%
  write_csv(path = "./OUTPUT/TABLES/yield_correlations_clump_spacing.csv") %>%
  dplyr::ungroup() %>%
  dplyr::summarise(mean(Correlation))

##### yield_correlations_clump_spacing.pdf #####
## A visual comparison between the scaled experimental and VEP estimated yields
## with linear regression
phi <- (1+sqrt(5))/2
mai <- c(0.25,0.25,0,0)
fig.width <- 5
fig.height <- fig.width/phi

quartz(file="./OUTPUT/FIGURES/yield_correlations_clump_spacing.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=6, dpi=600)
par(mai=mai, xpd=F)

yields %>%
  dplyr::select(Garden, Season, `Yield by clump area`) %>%
  dplyr::group_by(Garden, Season) %>%
  dplyr::summarise(`Yield by clump area` = mean(`Yield by clump area`)) %>%
  dplyr::full_join(PFP_VEP_yields, by = c("Garden","Season")) %>%
  dplyr::rename(`Experimental yield` = `Yield by clump area`, `Estimated yield` = Yield) %>%
  dplyr::mutate(`Experimental yield` = scale(`Experimental yield`) %>% as.numeric(), 
                `Estimated yield` = scale(`Estimated yield`) %>% as.numeric()) %>%
  ggplot(aes(x = `Estimated yield`, y = `Experimental yield`, color = Garden)) + 
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE,
              size = 0.5) +   # Don't add shaded confidence region
  geom_point() +
  scale_colour_brewer("Garden", palette="OrRd") +
  theme(axis.text=element_text(size=6, color = "black"),
        title=element_text(size=8,face="bold", color = "black"),
        legend.text=element_text(size=6, color = "black"))

dev.off()
distill("./OUTPUT/FIGURES/yield_correlations_clump_spacing.pdf")

##### END yield_correlations_clump_spacing.pdf #####

##### vepi_productivity.pdf #####
# A smoothing distribution to use throughout the analyses
# This is a 21-year wide Gaussian distribution with a 
# mean of 0 and standard deviation of 5 years
dist <- dnorm(seq(-10,10,1), sd=5)

# Import VEPI Paleoproductivity
VEPI.paleoprod <- raster::brick("./DATA/VEPI_paleoprod.tif")

VEPI.paleoprod.mean <- mean(VEPI.paleoprod[])
VEPI.paleoprod.mean.spatial <- mean(VEPI.paleoprod)
VEPI.paleoprod.mean.temporal <- cellStats(VEPI.paleoprod, mean)
VEPI.paleoprod.mean.temporal.smooth <- stats::filter(VEPI.paleoprod.mean.temporal,filter = dist)

plot.width <- 5
plot.height <- plot.width * (nrow(VEPI.paleoprod)/ncol(VEPI.paleoprod))

fig.width <- plot.width + 0.2
fig.height <- plot.height + 0.3 + 1.5

palette <- brewer.pal(11, "RdYlGn")
colors.begin <- rev(colorRampPalette(rev(palette[c(1:6)]),bias=1.2)(round(VEPI.paleoprod.mean)))
colors.mid <- colorRampPalette(palette[c(6:11)],bias=1.2)(600 - round(VEPI.paleoprod.mean))
colors.end <- rep(palette[11],1900)
colors <- c(colors.begin,colors.mid, colors.end)

quartz(file="./OUTPUT/FIGURES/vepi_productivity.pdf", width=fig.width, height=fig.height, antialias=FALSE, bg="white", type='pdf', family="Gulim", pointsize=8, dpi=600)
par(mai=c(1.5 + 0.2,0.1,0.1,0.1), xpd=F)

plot(1, type='n', xlab="", ylab="", 
     xlim=c(extent(VEPI.paleoprod)@xmin,extent(VEPI.paleoprod)@xmax),
     ylim=c(extent(VEPI.paleoprod)@ymin,extent(VEPI.paleoprod)@ymax), 
     xaxs="i", yaxs="i", axes=FALSE, main='')
plot(VEPI.paleoprod.mean.spatial, maxpixels=1000000, zlim=c(0,2500),add=T, col=colors, useRaster=TRUE, legend=FALSE)

par(mai=c(0.1,0.1,0.2 + plot.width * (nrow(VEPI.paleoprod)/ncol(VEPI.paleoprod)),0.1), xpd=T, new=T)
plot(1, type='n', xlab="", ylab="", xlim=c(0,5), ylim=c(0,600), xaxs="i", yaxs="i", axes=FALSE, main='')

legend.breaks <- seq(from=0, to=2500, length.out=(length(colors)+1))
rect(col=colors[1:600], border=NA, ybottom=0:599, ytop=1:600, xleft=0.15, xright=0.35, xpd=T)
# text(x = 0.5, y=0, labels=0, adj=c(0.5,0), cex=0.65, family='Helvetica Bold')
# text(x = 0.5, y=600, labels="> 600", adj=c(0.5,0.5), cex=0.65, family='Helvetica Bold')
# text(x = 0.5, y=round(VEPI.paleoprod.mean), labels=round(VEPI.paleoprod.mean), adj=c(0.5,0.5), cex=0.65, family='Helvetica Bold')
text(x = 0, y=300, labels="Yeild (kg/ha)", adj=c(0.5,1), cex=1, srt = 90, family='Helvetica Bold')
text(x = 0.5, y = c(0,100,200,round(VEPI.paleoprod.mean),300,400,500,600), labels = c(0,100,200,round(VEPI.paleoprod.mean),300,400,500,"> 600"), adj=c(0.5,0.5), cex=0.65, family='Helvetica Bold')

par(mai=c(0.1,0.85,0.2 + plot.width * (nrow(VEPI.paleoprod)/ncol(VEPI.paleoprod)),0.1), xpd=T, new=T)
plot(1, type='n', xlab="", ylab="", xlim=c(580,1320), ylim=c(0,600), xaxs="i", yaxs="i", axes=FALSE, main='')
abline(h = round(VEPI.paleoprod.mean),col = "gray50", lty = 1, xpd = F)
segments(x0 = seq(600,1300,100), x1 = seq(600,1300,100), y0 = 0, y1 = 500, col = "gray50", lty = 3)
lines(y = VEPI.paleoprod.mean.temporal, x = 600:1300)
lines(y = VEPI.paleoprod.mean.temporal.smooth, x = 600:1300, col = "red")
text(x = seq(600,1300,100), y = 550, labels = seq(600,1300,100), adj = c(0.5,1), cex=0.8, family='Helvetica Bold')
text(x = 950, y = 600, labels = "Year AD", adj = c(0.5,1), cex=1, family='Helvetica Bold')

axis(2, at = c(0,100,200,round(VEPI.paleoprod.mean),300,400,500,600), labels = F)
dev.off()
distill("./OUTPUT/FIGURES/vepi_productivity.pdf")
##### END vepi_productivity.pdf #####

##### Burns and Van West figures #####
file.copy("./DATA/BURNS.pdf","./OUTPUT/FIGURES/burns.pdf", overwrite = T)
distill("./OUTPUT/FIGURES/burns.pdf", gray = T)

file.copy("./DATA/VANWEST.pdf","./OUTPUT/FIGURES/vanwest.pdf", overwrite = T)
distill("./OUTPUT/FIGURES/vanwest.pdf", gray = T)
##### END Burns and Van West figures #####

##### Zip up SI data #####
zip("./Varien_Bocinsky_2017.zip", c("./DATA/","./Ethnobiology.R","./Ethnobiology.Rproj","./src/"), flags = "-r9X", extras = "-r --exclude=*.DS_Store* --exclude=*.git*")
##### END Zip up SI data #####

