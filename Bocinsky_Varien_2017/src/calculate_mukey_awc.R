calculate_mukey_awc <- function(SSURGO){
  # Horizon data
  # This is how the AWC is calculated!
  # This section recreates what is described in Johnson 2006:133--138
  # Calculating the AWC in the upper 6 and lower 54 inches of each soil
  
  chorizon <- dplyr::as_data_frame(SSURGO$tabular$chorizon) %>%
    dplyr::mutate(hzdept.r = hzdept.r * 0.393701, hzdepb.r = hzdepb.r * 0.393701) %>%  # Convert depth measurements to inches.
    dplyr::filter(hzdept.r < 60) %>% # Only include soil horizons that begin less than 60 inches below the surface
    dplyr::mutate(hzdepb.r = ifelse(hzdepb.r > 60,60,hzdepb.r)) %>% # Truncate horizons that end below 60 inches
    dplyr::filter(!is.na(awc.r)) %>% # Remove horizons with no AWC
    dplyr::group_by(cokey) # Group by cokey
  
  # If a horizon spans the 6-inch depth, split it in two at 6 inches 
  spans <- chorizon %>%
    dplyr::filter(hzdept.r<6 & hzdepb.r>6)
  spans <- dplyr::bind_rows(spans %>%
                              dplyr::mutate(hzdepb.r = ifelse(hzdepb.r > 6,6,hzdepb.r)), # Truncate horizons that end below 60 inches
                            spans %>%
                              dplyr::mutate(hzdept.r = ifelse(hzdept.r < 6,6,hzdept.r))
  ) %>%
    dplyr::arrange(cokey,chkey)
  chorizon <- dplyr::bind_rows(chorizon %>%
                                 dplyr::filter(hzdept.r>=6 | hzdepb.r<=6),
                               spans) %>%
    dplyr::arrange(cokey,hzdept.r) # Sort horizons by top depths and component keys
  
  
  # Join component AWCs to component table
  components <- dplyr::left_join(dplyr::as_data_frame(SSURGO$tabular$component),
                                chorizon %>%
                                  dplyr::mutate(awc = (hzdepb.r - hzdept.r) * awc.r, level = ifelse(hzdepb.r<=6,"upper_awc","lower_awc")) %>%  # Calculate final AWC for each horizon and Code horizons as being either 0-6 (level 1) or 6-60 (level 2) inches in depth
                                  dplyr::group_by(cokey,level) %>%
                                  dplyr::summarise(total_awc = sum(awc)) %>% # Sum AWS by level
                                  tidyr::spread(level,total_awc),
                                by = "cokey"
  )
  
  SSURGO$spatial@data <- dplyr::left_join(SSURGO$spatial@data,
                                          components %>%
                                            dplyr::group_by(mukey) %>%
                                            dplyr::summarise(upper_awc = sum(comppct.r * upper_awc, na.rm = T),lower_awc = sum(comppct.r * lower_awc, na.rm = T)) %>%
                                            dplyr::mutate(mukey = as.factor(mukey)),
                                          by = c("MUKEY" = "mukey")
  )
  
  return(SSURGO)
}
