library(tidyverse)

resample_transect <- function(x, randomseed = 1, nbootstrap = 10, nquad = 15, 
                              npoints = 25) {

  set.seed(randomseed)
  #First, we sum all the transects that belong to the same site.
  
  y <- x %>%
    group_by(Locality, Site, Loc, InsularCoast, Position) %>%
    summarise_if(is.numeric, sum) %>% 
    data.frame(stringsAsFactors = F, check.names = F)
  
  NewDF <- x[0, -1]
  NewDF$Transect <- as.character(NewDF$Transect)
  #For each row we will generate a repetition of each spp according to its frequency 
  
  for(i in 1:nrow(y)) {
    
    site_points <- NULL
    
    z <- y[i, y[i, ] != 0][, -c(1:5)] #exclusion of empty spp
  
    for(j in 1:length(z)) {
      w <- rep(names(z)[j], times = z[1, j])
      site_points <- c(site_points, w)
    }
  
    for(boot in 1:nbootstrap) {
      NewTransect <- x[0, -1]
      
      Factors <- y[i, 1:5]
      Factors$Transect <- paste("t", boot+4, sep = "") #we designate the new transects from 5 onwards
      Factors <- Factors[, c(1:2, 6, 3:5)] #We rearrange the factors order
      
      NewTransect[1, 1:6] <- Factors
      p <- sample(site_points, npoints*nquad, replace = F) %>% 
        table
      NewTransect[, names(p)] <- p
      NewDF <- rbind(NewDF, NewTransect)
    }
  }
  NewDF[is.na(NewDF)] <- 0
  return(NewDF)
}
