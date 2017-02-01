calculate_mukey_hand_planting_restriction <- function(SSURGO){
  ### HAND-PLANTING RESTRICTION BY MUKEY AND CELL ###
  ## This section takes the hand-planting classification from sdvattribute.txt 
  ## from the muaggatt tables and generates a raster of average restrictions per simulation cell.
  # Load USDA NRCS data.
  # Load the USGS survey data from the map units specified.
  # Load the hand-planting data
  SCM_RED <- dplyr::as_data_frame(SSURGO$tabular$cointerp) %>%
    dplyr::filter(mrulename=='FOR - Hand Planting Suitability') %>%
    dplyr::select(cokey,interpll) %>%
    dplyr::mutate(interpll = ifelse(is.na(interpll),0,interpll)) %>% # Set NAs to 0.0 reduction
    dplyr::group_by(cokey) %>%
    dplyr::summarise(SCM_RED = 1-(max(interpll, na.rm = T)/2)) # Take the maximum reduction per component
    # The reduction multiplier is 1 minus half the reduction value.

  # Join component SCM_REDs to component table
  components <- dplyr::left_join(dplyr::as_data_frame(SSURGO$tabular$component),
                                 SCM_RED,
                                 by = "cokey"
  )
  
  SSURGO$spatial@data <- dplyr::left_join(SSURGO$spatial@data,
                                          components %>%
                                            dplyr::group_by(mukey) %>%
                                            dplyr::summarise(SCM_RED = sum(comppct.r * SCM_RED, na.rm = T)) %>%
                                            dplyr::mutate(mukey = as.factor(mukey)),
                                          by = c("MUKEY" = "mukey")
  )
  
  return(SSURGO)
}
