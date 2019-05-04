BootResampler <- function(DataFrame, nbTransects = 10, nbQuadrats = 15, nbPoints = 25,
                          transects_folder = "data/transects/", SpMatrix_Path = FALSE,
                          randomseed = 1) {
  if(!is.data.frame(SpMatrix_Path)) {
    stop("Debes agregar la matriz de especies al parametro SpMatrix_Path, como data.frame")
  }
  source('r/RandomChoose.R')
  transect_files = list.files(transects_folder)
  Especies <- SpMatrix_Path
  
  BootDF <- DataFrame[0, ]
  l = 1
  set.seed(randomseed)
  
  for(i in 1:length(unique(DataFrame$Site))) {
    #if(l >= 3) next()
    y <- transect_files[gsub("_.*", "", transect_files) == unique(gsub("_.*", "", DataFrame$File))[i]] #Número de transectas en el sitio i
    
    Factors <- DataFrame[grep(unique(gsub("_.*", "", DataFrame$File))[i], DataFrame$File), 2:7][1, ] #Extracción de los factores correspondientes al sitio
    
    bootTransects <- rep(y, ceiling(nbTransects/length(y)))[1:nbTransects] # Repite las transectas disponibles de forma sucesiva, de manera que se alcance el N de transectas deseado
    
    for(j in 1:nbTransects) {
      #if(l >= 3) next()
      transect <- read.csv(paste("data/transects/", bootTransects[j], sep = ""), stringsAsFactors = F)[, 1:9]
      
      for(k in 1:length(unique(transect$spp.ID))) {
        transect[transect$spp.ID == unique(transect$spp.ID)[k], 3] <- as.character(Especies[, 1][Especies$Species.ID %in% unique(transect$spp.ID)[k]])
      }
      
      transect_boot <- randomchooseNpoints(transect, nquadrats = nbQuadrats, npoints = nbPoints)
      
      vec_sum <- aggregate(N.pts.per.species ~ spp.Name, data = transect_boot, FUN = sum)
      
      if(length(which(!vec_sum$spp.Name %in% names(BootDF))) > 0) { #Agrega columnas si hacen falte en BootDF, dado que no estén en la matriz original
        NewCols <- vec_sum$spp.Name[which(!vec_sum$spp.Name %in% names(BootDF))]
        BootDF <- cbind(BootDF, array(dim = c(nrow(BootDF), length(NewCols)), dimnames = list(NULL, NewCols)))
      }
      
      logi_transect <- which(names(BootDF) %in% vec_sum$spp.Name)
      logi_transect <- logi_transect[order(names(BootDF)[logi_transect], decreasing = F)]
      
      BootDF[l, logi_transect] <- vec_sum$N.pts.per.species
      
      BootDF[l, 1] <- bootTransects[j]
      BootDF[l, 2:7] <- Factors
      BootDF[l, 4] <- paste("t_boot_", j, sep = "")
      BootDF[is.na(BootDF)] <- 0
      
      l <- l + 1
      
    }
    
  }
  
  BootDF
}
