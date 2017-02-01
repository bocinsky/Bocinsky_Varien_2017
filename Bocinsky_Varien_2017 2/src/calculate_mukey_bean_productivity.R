calculate_mukey_bean_productivity <- function(SSURGO){
  beans <- dplyr::as_data_frame(SSURGO$tabular$cocropyld) %>%
    dplyr::filter(grepl("bean",cropname)) %>%
    dplyr::select(cokey,nonirryield.r) %>%
    dplyr::rename(bean_yield = nonirryield.r)
    
  components <- dplyr::left_join(dplyr::as_data_frame(SSURGO$tabular$component),
                                 beans,
                                 by = "cokey"
  )
  
  SSURGO$spatial@data <- dplyr::left_join(SSURGO$spatial@data,
                                          components %>%
                                            dplyr::group_by(mukey) %>%
                                            dplyr::summarise(bean_yield = sum(comppct.r * bean_yield, na.rm = T)) %>%
                                            dplyr::mutate(mukey = as.factor(mukey)),
                                          by = c("MUKEY" = "mukey")
  )
  
  return(SSURGO)
  
}
