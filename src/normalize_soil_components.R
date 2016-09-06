normalize_soil_components <- function(SSURGO){
  components <- dplyr::as_data_frame(SSURGO$tabular$component) %>%
    dplyr::group_by(mukey)
  
  components <- components %>%
    dplyr::left_join(
      components %>%
        dplyr::summarise(comppct.r.sum = sum(comppct.r)
        ),
      by = "mukey") %>%
    dplyr::mutate(comppct.r = comppct.r/comppct.r.sum) %>%
    dplyr::select(-comppct.r.sum)
  
  if("data.frame" %in% class(SSURGO$tabular$component) & !("tbl_df" %in% class(SSURGO$tabular$component))){
    SSURGO$tabular$component <- as.data.frame(components)
  }else{
    SSURGO$tabular$component <- components
  }
  
  return(SSURGO)
}
