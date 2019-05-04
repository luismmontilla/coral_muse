library(tidyverse)

resample_quadrat <- function(x, randomseed = 1, nbootstrap = 10, nquads = 15) {
  
  set.seed(randomseed)
  #First, we sum all the transects that belong to the same site.
  
  y <- x %>%
    group_by(Locality, Site, Loc, InsularCoast, Position, Transect, Quadrat) %>%
    summarise_if(is.numeric, sum) %>% 
    data.frame(stringsAsFactors = F, check.names = F)
  
  NewDF <- x[0, -1]
  NewDF$Transect <- as.character(NewDF$Transect)
  NewDF <- NewDF[, -4]
  #For each row we will generate a repetition of each spp according to its frequency 
  
  for(i in 1:length(unique(y$Site))) {
    
    site_quads <- y[y$Site == unique(y$Site)[i], ]
    
    for(boot in 1:nbootstrap) {
      NewTransect <- NewDF[0, ]
      
      Factors <- site_quads[1, 1:5]
      Factors$Transect <- paste("t", boot+4, sep = "") #we designate the new transects from 5 onwards
      Factors <- Factors[, c(1:2, 6, 5:3)] #We rearrange the factors order
      
      NewTransect[1, 1:6] <- Factors
      p <- site_quads[sample(1:nrow(site_quads), nquads, replace = T), ]
      p <- colSums(p[, c(8:ncol(p))])
      
      NewTransect[, names(p)] <- p
      NewDF <- rbind(NewDF, NewTransect)
    }
  }
  NewDF[is.na(NewDF)] <- 0
  return(NewDF)
}
