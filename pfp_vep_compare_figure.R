library(raster)
library(sf)
library(magrittr)
library(villager)

gardens <-
  sf::read_sf("https://raw.githubusercontent.com/crowcanyon/pfp_ebook/master/data/gardens.geojson")

yields <-
  readr::read_csv("https://pfp.crowcanyon.org/docs/tables/yield_table.csv")

# daymet <-
#   gardens %>%
#   FedData::get_daymet(template = .,
#                       label = "PFP Gardens",
#                       elements = c("prcp", "tmax", "tmin"),
#                       years = 1981:2019,
#                       tempo = "mon",
#                       extraction.dir = "DATA/daymet/"
#   )
# 
# gardens_daymet <-
#   daymet %>%
#   purrr::map_dfr(~raster::extract(.x, 
#                                   y = gardens %>%
#                                     sf::st_centroid(),
#                                   df = TRUE) %>%
#                    dplyr::mutate(ID = gardens$Name) %>%
#                    dplyr::rename(Name = ID),
#                  .id = "element") %>%
#   tidyr::pivot_longer(-element:-Name, names_to = "date") %>%
#   dplyr::mutate(date = date %>%
#                   stringr::str_remove("X") %>%
#                   lubridate::as_date(),
#                 Year = lubridate::year(date),
#                 Month = lubridate::month(date)) %>%
#   tidyr::pivot_wider(names_from = element,
#                      values_from = value) %>%
#   dplyr::mutate(prcp = units::as_units(prcp,"mm") %>% units::set_units("in"),
#                 tmax = units::as_units(tmax,"deg_c") %>% units::set_units("deg_f"),
#                 tmin = units::as_units(tmin,"deg_c") %>% units::set_units("deg_f"),
#                 tavg = (tmax + tmax) / 2 ) %>%
#   dplyr::select(-tmax, -tmin, -date) %>%
#   tidyr::nest(daymet = c(-Name)) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(normals = 
#                   list(
#                     daymet %>%
#                       dplyr::filter(Year %in% 1981:2010) %>%
#                       dplyr::group_by(Month) %>%
#                       dplyr::summarise_at(.vars = dplyr::vars(prcp, tavg), mean)
#                   )
#   )

gardens_ssurgo <- 
  gardens %>%
  sf::st_centroid() %>% 
  sf::st_buffer(dist = units::as_units(0.00001,"°")) %>%
  dplyr::rowwise() %>%
  dplyr::group_split() %>%
  magrittr::set_names(gardens$Name) %>%
  purrr::map(sf::st_as_sf) %>%
  purrr::map(function(x){
    FedData::get_ssurgo(template = x %>%
                          as("Spatial"),
                        label = x$Abbreviation,
                        extraction.dir = "DATA/ssurgo/",
                        raw.dir = "DATA/ssurgo/raw/"
    )
  }) %>%
  purrr::map_dfr(function(x){
    villager:::calculate_mukey_awc(x) %$%
      spatial %>%
      sf::st_as_sf() %>%
      dplyr::select(lower_awc) %>%
      sf::st_drop_geometry() %>%
      dplyr::mutate(lower_awc = lower_awc / 100)
  },
  .id = "Name")

# Use FedData to download the weather data from the Cortez GHCN
cortez_prcp <- 
  FedData::get_ghcn_daily_station(ID="USC00051886", 
                                  elements = c("PRCP"), 
                                  raw.dir = "./OUTPUT/DATA/GHCN") %>%
  FedData::station_to_data_frame() %>% # convert station data matrix to a data frame
  tibble::as_tibble() %>%
  dplyr::mutate(PRCP = units::set_units(PRCP / 10, "mm") %>% units::set_units("in"),
                Season = lubridate::year(DATE),
                Season = ifelse(lubridate::month(DATE) %in% 10:12, Season + 1, Season)) %>%
  dplyr::group_by(Season) %>%
  dplyr::summarise(`Precipitation (in.)` = sum(PRCP, na.rm = T))



library(ggplot2)
plot_data <- 
  yields %>%
  dplyr::left_join(gardens %>%
                     dplyr::select(Name, Abbreviation) %>%
                     sf::st_drop_geometry(),
                   by = c("Garden"="Abbreviation")) %>%
  dplyr::left_join(cortez_prcp) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Season = 
                  paste0(Season,"\n",round(`Precipitation (in.)`,2),'"') %>%
                  factor(ordered = TRUE),
                Name = factor(Name, ordered = TRUE),
                Garden = factor(Garden, 
                                levels = c("PHG", 
                                           "KUG", 
                                           "PLC", 
                                           "POG", 
                                           "CDG", 
                                           "MCG"),
                                ordered = TRUE)
  )

plot_data %>%
  ggplot(aes(x = Season)) +
  geom_line(aes(x = Season,
                y = units::drop_units(`Precipitation (in.)`) * 200,
                group = 1),
            data = 
              plot_data %>%
              dplyr::select(Season, `Precipitation (in.)`) %>%
              dplyr::distinct(),
            color = "black",
            linetype = 1,
            size = 1) +
  geom_boxplot(aes(y = `PFP experimental yield (kg/ha)`,
                   fill = Name),
               outlier.shape = NA,
               varwidth = FALSE,
               position = position_dodge2(padding = 0.3,
                                         preserve = "single")) +
  # geom_line(size = 2) +
  
  scale_fill_brewer(type = "qual",
                    palette = "Set2") +
  scale_y_continuous(breaks = seq(0,4000,500),
                     limits = c(0,4000),
                     expand = expansion(0,0),
                     sec.axis = sec_axis(
                       trans = ~ . / 200,
                       name = "Oct–Sept Precipitation (in.)"
                     )) +
  xlab(NULL) +
  ylab("PFP yield (kg/ha)") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("OUTPUT/FIGURES/pfp_yields_box.pdf",
       width = 10.5,
       height = 6,
       device = cairo_pdf)


plot_data %>%
  ggplot(aes(x = Season)) +
  geom_line(aes(x = Season,
                y = units::drop_units(`Precipitation (in.)`) * 200,
                group = 1),
            data = 
              plot_data %>%
              dplyr::select(Season, `Precipitation (in.)`) %>%
              dplyr::distinct(),
            color = "black",
            linetype = 1,
            size = 1) +
  geom_boxplot(aes(y = `PFP experimental yield (kg/ha)`,
                   fill = Garden),
               outlier.shape = NA,
               varwidth = FALSE,
               position = position_dodge2(padding = 0.3,
                                          preserve = "single")) +
  # geom_line(size = 2) +
  
  scale_fill_brewer(type = "qual",
                    palette = "Set2") +
  scale_y_continuous(breaks = seq(0,4000,500),
                     limits = c(0,4000),
                     expand = expansion(0,0),
                     sec.axis = sec_axis(
                       trans = ~ . / 200,
                       name = "Oct–Sept Precipitation (in.)"
                     )) +
  xlab(NULL) +
  ylab("PFP yield (kg/ha)") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("OUTPUT/FIGURES/pfp_yields_box_side.pdf",
       width = 10.5,
       height = 6,
       device = cairo_pdf)

ggsave("OUTPUT/FIGURES/pfp_yields_box_side.png",
       width = 10.5,
       height = 6)



plot_data <- 
  yields %>%
  dplyr::group_by(Season,Garden) %>%
  dplyr::summarise(`PFP experimental yield (kg/ha)` = mean(`PFP experimental yield (kg/ha)`)) %>%
  dplyr::left_join(gardens %>%
                     dplyr::select(Name, Abbreviation) %>%
                     sf::st_drop_geometry(),
                   by = c("Garden"="Abbreviation")) %>%
  dplyr::left_join(cortez_prcp) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Season = 
                  paste0(Season,"\n",round(`Precipitation (in.)`,2),'"') %>%
                  factor(ordered = TRUE),
                Name = factor(Name, ordered = TRUE),
                Garden = factor(Garden, ordered = TRUE)
  )

plot_data %>%
  ggplot(aes(x = Season)) +
  geom_line(aes(x = Season,
                y = units::drop_units(`Precipitation (in.)`) * 100,
                group = 1),
            data = 
              plot_data %>%
              dplyr::select(Season, `Precipitation (in.)`) %>%
              dplyr::distinct(),
            color = "black",
            linetype = 2,
            size = 1) +
  # geom_boxplot(aes(y = `PFP experimental yield (kg/ha)`,
  #                  fill = Name),
  #              outlier.shape = NA) +
  geom_line(aes(y = `PFP experimental yield (kg/ha)`,
                color = Name,
                group = Name),
            size = 1.5) +
  scale_color_brewer(type = "qual",
                    palette = "Set2") +
  scale_y_continuous(breaks = seq(0,3000,500),
                     limits = c(0,3000),
                     expand = expansion(0,0),
                     sec.axis = sec_axis(
                       trans = ~ . / 100,
                       name = "Oct–Sept Precipitation (in.)"
                     )) +
  xlab(NULL) +
  ylab("PFP yield (kg/ha)") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("OUTPUT/FIGURES/pfp_yields_lines.pdf",
       width = 10.5,
       height = 6,
       device = cairo_pdf)


plot_data %>%
  ggplot(aes(x = Season)) +
  geom_line(aes(x = Season,
                y = units::drop_units(`Precipitation (in.)`) * 100,
                group = 1),
            data = 
              plot_data %>%
              dplyr::select(Season, `Precipitation (in.)`) %>%
              dplyr::distinct(),
            color = "black",
            linetype = 2,
            size = 1) +
  # geom_boxplot(aes(y = `PFP experimental yield (kg/ha)`,
  #                  fill = Name),
  #              outlier.shape = NA) +
  geom_line(aes(y = `PFP experimental yield (kg/ha)`,
                color = Garden,
                group = Garden),
            size = 1.5) +
  scale_color_brewer(type = "qual",
                     palette = "Set2") +
  scale_y_continuous(breaks = seq(0,3000,500),
                     limits = c(0,3000),
                     expand = expansion(0,0),
                     sec.axis = sec_axis(
                       trans = ~ . / 100,
                       name = "Oct–Sept Precipitation (in.)"
                     )) +
  xlab(NULL) +
  ylab("PFP yield (kg/ha)") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("OUTPUT/FIGURES/pfp_yields_lines_side.pdf",
       width = 10.5,
       height = 6,
       device = cairo_pdf)



plot_data %>%
  ggplot(aes(x = Season)) +
  geom_line(aes(x = Season,
                y = units::drop_units(`Precipitation (in.)`) * 100,
                group = 1),
            data = 
              plot_data %>%
              dplyr::select(Season, `Precipitation (in.)`) %>%
              dplyr::distinct(),
            color = "black",
            linetype = 2,
            size = 1) +
  # geom_boxplot(aes(y = `PFP experimental yield (kg/ha)`,
  #                  fill = Name),
  #              outlier.shape = NA) +
  geom_point(aes(y = `PFP experimental yield (kg/ha)`,
                color = Name,
                group = Name),
            size = 1.5) +
  scale_color_brewer(type = "qual",
                    palette = "Set2") +
  scale_y_continuous(breaks = seq(0,3000,500),
                     limits = c(0,3000),
                     expand = expansion(0,0),
                     sec.axis = sec_axis(
                       trans = ~ . / 100,
                       name = "Oct–Sept Precipitation (in.)"
                     )) +
  xlab(NULL) +
  ylab("PFP yield (kg/ha)") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("OUTPUT/FIGURES/pfp_yields_points.pdf",
       width = 10.5,
       height = 6,
       device = cairo_pdf)


plot_data %>%
  ggplot(aes(x = Season)) +
  geom_line(aes(x = Season,
                y = units::drop_units(`Precipitation (in.)`) * 100,
                group = 1),
            data = 
              plot_data %>%
              dplyr::select(Season, `Precipitation (in.)`) %>%
              dplyr::distinct(),
            color = "black",
            linetype = 2,
            size = 1) +
  # geom_boxplot(aes(y = `PFP experimental yield (kg/ha)`,
  #                  fill = Name),
  #              outlier.shape = NA) +
  geom_point(aes(y = `PFP experimental yield (kg/ha)`,
                 color = Garden,
                 group = Garden),
             size = 1.5) +
  scale_color_brewer(type = "qual",
                     palette = "Set2") +
  scale_y_continuous(breaks = seq(0,3000,500),
                     limits = c(0,3000),
                     expand = expansion(0,0),
                     sec.axis = sec_axis(
                       trans = ~ . / 100,
                       name = "Oct–Sept Precipitation (in.)"
                     )) +
  xlab(NULL) +
  ylab("PFP yield (kg/ha)") +
  theme_minimal(18) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("OUTPUT/FIGURES/pfp_yields_points_side.pdf",
       width = 10.5,
       height = 6,
       device = cairo_pdf)


# test <- 
#   gardens %>%
#   dplyr::left_join(gardens_daymet) %>%
#   dplyr::left_join(gardens_ssurgo) %>%
#   dplyr::rowwise() %>%
#   dplyr::mutate(pdsi = 
#                   list(list(
#                     villager::scpdsi(
#                       monthly_T = daymet$tavg, 
#                       monthly_P = daymet$prcp, 
#                       mon_T_normal = normals$tavg, 
#                       awc = lower_awc, 
#                       lat = geometry %>%
#                         sf::st_centroid() %>%
#                         sf::st_coordinates() %>%
#                         magrittr::extract(,"Y"), 
#                       start_year = min(daymet$Year),                         
#                       end_year = max(daymet$Year)
#                     )
#                 ))
#   )
# 
# villager::scpdsi(monthly_T = 70:92, 
#                  monthly_P = 1:12, 
#                  mon_T_normal = 70:92, 
#                  awc = 7, 
#                  lat = 45, 
#                  start_year = 1990, 
#                  end_year = 1990)
# 
# 
# 
# 
# gardens %$%
#   geometry %>%
#   sf::st_centroid() %>%
#   sf::st_coordinates() %>%
#   magrittr::extract(,"Y")
# 
# 
# ## Produce PDSI estimates for each soil
# # Download the scpdsi C code from GitHub
# FedData::download_data("https://raw.githubusercontent.com/village-ecodynamics/scpdsi/master/src/scpdsi.cpp",
#                        destdir = "./src/")
# 
# # Load the PDSI c program
# Rcpp::sourceCpp("./src/scpdsi.cpp")
# 
# 
# 
# # Run the PDSI program for each soil
# # This operates in parallel; it will register the number of cores appropriate to your machine
# registerDoParallel()
# vepi.soils.ccac.pdsi <- foreach::foreach(soil = 1:nrow(vepi.soils.ccac)) %dopar% {
#   
#   out <- scpdsi(monthly_T = cortez_weather_monthly$TAVG_F,
#                 monthly_P = cortez_weather_monthly$PRCP_IN,
#                 mon_T_normal = cortez_weather_monthly_norms$TAVG_F,
#                 awc = vepi.soils.ccac@data[soil,"AWC_Lower_median"],
#                 lat = vepi_latitude,
#                 start_year = burnin_start,
#                 end_year = max(seasons))
#   expand.grid(YEAR = burnin_start:max(seasons), MONTH = 1:12) %>%
#     tibble::as_tibble() %>%
#     dplyr::arrange(YEAR,MONTH) %>%
#     dplyr::mutate(PDSI = out) %>%
#     dplyr::filter(YEAR %in% seasons, MONTH == 6) %>%
#     dplyr::select(YEAR,PDSI) 
#   
# }
# stopImplicitCluster()
# 
# 
