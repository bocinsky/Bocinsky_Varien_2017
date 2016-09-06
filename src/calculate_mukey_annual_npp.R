calculate_mukey_annual_npp <- function(SSURGO){
  SSURGO <- normalize_soil_components(SSURGO)
  
  SSURGO$spatial@data <- dplyr::left_join(SSURGO$spatial@data,
                   dplyr::as_data_frame(SSURGO$tabular$component) %>%
                     dplyr::group_by(mukey) %>%
                     dplyr::summarise(NPP = sum(comppct.r * rsprod.r, na.rm = T)) %>%
                     dplyr::mutate(mukey = as.factor(mukey)),
                   by = c("MUKEY" = "mukey")
                   )

  return(SSURGO)
}