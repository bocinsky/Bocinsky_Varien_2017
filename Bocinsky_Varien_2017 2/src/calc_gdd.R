calc_gdd <- function(tmin, tmax, t.base, t.cap=NULL, multiplier=1, to_fahrenheit=T){
  if(length(tmin)!=length(tmax)){
    stop("tmin and tmax  must have same length!")
  }
  
  t.base <- t.base*multiplier
  if(!is.null(t.cap)){
    t.cap <- t.cap*multiplier
  }
    
    # Floor tmax and tmin at Tbase
    tmin[tmin < t.base] <- t.base
    tmax[tmax < t.base] <- t.base
    
    # Cap tmax and tmin at Tut
    if(!is.null(t.cap)){
      tmin[tmin > t.cap] <- t.cap
      tmax[tmax > t.cap] <- t.cap
    }
    
    GDD <- ((tmin+tmax)/2)-t.base

    if(to_fahrenheit){
      GDD <- GDD * 1.8
    }
  
  return(GDD)
}