library(FedData)
FedData::pkg_test("geojsonio")
FedData::pkg_test("readr")
FedData::pkg_test("raster")
FedData::pkg_test("readxl")

## Download the VEP I simulation data
FedData::download_data("https://github.com/crowcanyon/vep_sim_beyondhooperville/raw/master/VEPI_data.zip",
                       destdir = "./DATA/")
scmred <- unz("./DATA/VEPI_data.zip",
    filename = "VEPI_data/SCMRed.data") %>%
  readLines() %>%
  as.numeric() %>%
  tibble::as_tibble() %>%
  dplyr::mutate(VPSCode = row_number()) %>%
  dplyr::rename(SCM_RED = value)

unz("./DATA/VEPI_data.zip",
    filename = "VEPI_data/newsoil.data") %>%
  scan()

# Load the VEPI soils information
vepi <- readr::read_csv("./DATA/village_soil.csv") %>%
  dplyr::left_join(scmred, by = "VPSCode")

# Load a polygon of the Indian Camp Ranch + Crow Canyon study area
CCAC <- rgdal::readOGR(dsn = "./DATA/ccac.geojson", "OGRGeoJSON", verbose = FALSE)

# Download the NRCS SSURGO soils data for the study area
# This downloads the soils data for the CO671 soil survey
vepi_soils_ccac <- FedData::get_ssurgo(template = CCAC,
                                       label = "CCAC",
                                       raw.dir = "./OUTPUT/DATA/SSURGO/RAW/",
                                       extraction.dir = "./OUTPUT/DATA/SSURGO/EXTRACTIONS/",
                                       force.redo =  F)
# Add the data to the shapefile of soil mapunits
vepi_soils_ccac$spatial@data <- dplyr::left_join(vepi_soils_ccac$spatial@data,vepi_soils_ccac$tabular$mapunit %>%
                                 dplyr::select(mukey,muname) %>%
                                 dplyr::mutate(mukey = as.factor(mukey)),
                               by = c("MUKEY" = "mukey"))
vepi_soils_ccac <- vepi_soils_ccac$spatial


# Join the VEPI soils data to the SSURGO spatial dataset
vepi_soils_ccac@data %<>%
  dplyr::select(MUSYM,muname,MUKEY) %>%
  dplyr::mutate(MUSYM = as.character(MUSYM)) %>%
  dplyr::left_join(vepi, by = c("MUSYM" = "Musym")) %>%
  dplyr::mutate(MUSYM = as.numeric(MUSYM)) %>%
  dplyr::select(MUSYM,MUKEY,VPSCode,CLUSTER,SCM_RED,muname) %>%
  as.data.frame() %>%
  dplyr::left_join(vepi %>%
                     dplyr::group_by(CLUSTER) %>%
                     dplyr::summarise(AWC_Lower_median = median(AWC_Lower, na.rm = T),
                                      NYProd_lb_ac = sum(X__of_Cluster * NYProd_lb_ac, na.rm = T)), by = "CLUSTER")

geojsonio::geojson_write(vepi_soils_ccac,file = "./DATA/vepi_soils_ccac.geojson")

vepi_soils_ccac <- rgdal::readOGR(dsn = "./DATA/vepi_soils_ccac.geojson", "OGRGeoJSON", verbose = FALSE)

spplot(vepi_soils_ccac, zcol = "CLUSTER")
spplot(vepi_soils_ccac, zcol = "AWC_Lower_median")
spplot(vepi_soils_ccac, zcol = "NYProd_lb_ac")
spplot(vepi_soils_ccac, zcol = "SCM_RED")


### Generate VEPI paleoprod raster
raster::rasterOptions(chunksize=2e+08,maxmemory=2e+09)

# Load a polygon of the VEPI study area
vepi_polygon <- FedData::polygon_from_extent(raster::extent(676400, 721800, 4126200, 4166200), 
                                             "+proj=utm +datum=NAD27 +zone=12")

# The VEP project uses a 200m grid. Set up the grid as a raster object.
vepi_raster <- raster::raster(ext = raster::extent(vepi_polygon),
                              nrows = abs(ymax(vepi_polygon)-ymin(vepi_polygon))/200,
                              ncols = abs(xmax(vepi_polygon)-xmin(vepi_polygon))/200,
                              crs = sp::CRS(raster::projection(vepi_polygon)))

# Import VEPI Paleoprod
VEPI.paleoprod.raw <- unz("./DATA/VEPI_data.zip",
    filename = "VEPI_data/al_year600-1300.dat") %>%
  scan()

# Create a raster brick of scaled MVDF series
VEPI.paleoprod.raw <- raster::brick(vepi_raster,
                                    values = F,
                                    nl = length(600:1300)) %>%
  raster::setValues(VEPI.paleoprod.raw)

# Inport correction values
# From Cell.java, line 2338:
# setColdCorrelation();//JAC 9/25/04 calculates cold correlation
# maize_pot = (int) (((veg * 10) + 4) * 2.36775 * adjust_factor * cold_corr * scmr);//JAC 9/25/04 adjusts maize production based on production file and cold correction
# where veg is the raw VEPI paleoprod reconstruction; adjust_factor = 0.682; cold_corr is a raster as a function of elevation; and scmr is the hand planting reduction raster
# set the adjust factor
adjust_factor <- 0.682

# Generate the cold correction raster
vepi.dem <- unz("./DATA/VEPI_data.zip",
                filename = "VEPI_data/dem.data") %>%
  scan()
VEPI.dem.rast <- vepi_raster %>%
  raster::setValues(vepi.dem)

# Load the Almagre and Mesa Verde Douglas Fir tree ring series used in VEP I (Greybill 1983; Dean and Robinson 1978:29–30)
FedData::download_data("https://raw.githubusercontent.com/crowcanyon/vep_paleoprod/master/DATA/vepi_tree_rings.csv",
                       destdir = "./DATA/")
vepi_tree_rings <- readr::read_csv("./DATA/vepi_tree_rings.csv", col_types = "idd") %>%
  dplyr::mutate(`Almagre Mountain B index scaled` = scale(`Almagre Mountain B index`) %>% as.vector(),
                `Mesa Verde Douglas Fir index scaled` = scale(`Mesa Verde Douglas Fir index`) %>% as.vector()) %>%
  dplyr::filter(Year %in% 600:1300)

coldcorr_upper_limit <- 2395 # Above which there is no production
coldcorr_lower_limit <- 2150 # Below which production isn't affected
min_almagre_z <- min(vepi_tree_rings$`Almagre Mountain B index scaled`, na.rm = T)

# Create a raster brick of scaled Almagre series
vepi_raster_Almagre_Z <- raster::brick(vepi_raster,
                                       values = F,
                                       nl = length(vepi_tree_rings$`Almagre Mountain B index scaled`)) %>%
  raster::setValues(1) %>%
  magrittr::multiply_by(vepi_tree_rings$`Almagre Mountain B index scaled`)

# Apply the cold-correction scaling function
vepi_raster_coldcorr <- {(coldcorr_upper_limit - VEPI.dem.rast)/(coldcorr_upper_limit-coldcorr_lower_limit)} * 
{(abs(min_almagre_z) + vepi_raster_Almagre_Z)/abs(min_almagre_z)}

# No change when Almagre is above the mean
vepi_raster_coldcorr[vepi_raster_Almagre_Z >=0] <- 1

rm(vepi_raster_Almagre_Z)
gc();gc()


# Disallow production over 2395 meters
vepi_raster_coldcorr[VEPI.dem.rast>=coldcorr_upper_limit] <- 0
# No changes under 2150 meters
vepi_raster_coldcorr[VEPI.dem.rast<=coldcorr_lower_limit] <- 1

# Don't allow decreases in production below zero, or increases in production!
vepi_raster_coldcorr[vepi_raster_coldcorr<0] <- 0
vepi_raster_coldcorr[vepi_raster_coldcorr>1] <- 1

# get handplanting reduction raster
soils <- unz("./DATA/VEPI_data.zip",
             filename = "VEPI_data/newsoil.data") %>%
  scan()
VEPI.soils.rast <- vepi_raster %>%
  raster::setValues(soils)
names(VEPI.soils.rast) <- "VPSCode"
CLUSTER.rast <- raster::subs(VEPI.soils.rast, vepi %>% dplyr::select(VPSCode,CLUSTER), by = "VPSCode", which = "CLUSTER")

# Generate the handplant reduction raster
SCMR.rast <- raster::subs(VEPI.soils.rast, vepi %>% dplyr::select(VPSCode,SCM_RED), by = "VPSCode", which = "SCM_RED")
SCMR.rast[is.na(SCMR.rast)] <- 0

#apply the potential productivity corrections
VEPI.paleoprod <- ((VEPI.paleoprod.raw * 10) + 4) * 2.36775 * adjust_factor * vepi_raster_coldcorr * SCMR.rast

raster::writeRaster(VEPI.paleoprod, "./DATA/vepi_paleoprod.tif",
                    datatype = "FLT4S",
                    options = c("COMPRESS=DEFLATE", "ZLEVEL=9", "INTERLEAVE=BAND"),
                    overwrite = T,
                    setStatistics = FALSE)

unlink("./DATA/VEPI_data.zip")
unlink("./DATA/VEPI_data/", recursive = T)


# sas7bdat::read.sas7bdat("~/git/vep_paleoprod/DATA/SAS/village_soil.sas7bdat") %>%
#   tibble::as_tibble() %>%
#   dplyr::select(Musym, VPSCode, CLUSTER, AWC_Lower, NYProd_lb_ac, X__of_Cluster,HECTARES) %>%
#   dplyr::filter(!is.nan(CLUSTER)) %>%
#   dplyr::arrange(Musym, -HECTARES) %>%
#   dplyr::distinct(Musym, .keep_all = TRUE) %>%
#   dplyr::select(-HECTARES) %>%
#   dplyr::arrange(VPSCode) %>%
#   readr::write_csv("./DATA/village_soil.csv")
