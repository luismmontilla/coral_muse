ErrorCheck <- data.frame(array(dim = c(0, 10)))

names(ErrorCheck) <- c("File", "Locality", "Site", "Transect", "Observer", "Quadrats",  "Total_points", "Unassigned", "ExtraAppend_Quadrats", "Unassigned_Quadrats")

#Lectura de las matrices correspondientes a cada transecta
Transectas <- list.files(path = 'data/transects/', 
                         pattern = "rndpts.csv",
                         full.names = TRUE) #vector que agrupa todos los nombres de archivos tipo rndots.csv que hay en la carpeta

Localidades <- read.csv("data/Venezuela_localities.csv", 
                        header = T, stringsAsFactors = F)

for(i in 1:length(Transectas)) {
  
  Transecta <- read.csv(Transectas[i], stringsAsFactors = F)
  Transecta <- Transecta[, 1:9]
  Nombres_col <- c("Image", "Path", "spp.Name", "spp.ID", "N.pts.per.species", "Cov..per.species", "N.pts.ALL.species", "Cov..ALL.species", "Reference.random.pts")
  
  if(all(names(Transecta) != Nombres_col)) {
    Transecta <- read.csv(Transectas[i], stringsAsFactors = F, header = F)
    Transecta <- Transecta[, 1:9]
    names(Transecta) <- Nombres_col
    
  }
  Transecta$spp.Name[grep("Ind", Transecta$spp.Name)] <- "Biotic (unknown)"
  
  #Extraer factores
  
  #Factores <- c(Transectas[i], "", strsplit((gsub(".rndpts.csv", "", Transectas[i])), split = "_")[[1]])
  Factores <- c(Transectas[i], "", strsplit((gsub(".rndpts.csv", "", Transectas[i])), split = "_")[[1]])
  Factores[2:3] <- as.character(Localidades[grep(gsub(".rndpts.csv", "", Transectas[i]), Localidades$File), ])[c(2:3)]
  
  
  #### Matriz para Verificar errores en numero de puntos ####
  if(length(Factores) == 5) {
    Factores[4] <- Factores[5]
  }
  ErrorCheck[i, 1:4] <- Factores
  ErrorCheck[i, 5]   <- as.character(Localidades[grep(gsub(".rndpts.csv", "", Transectas[i]), Localidades$File), ])[c(5)]
  ErrorCheck[i, 6]   <- length(unique(Transecta$Image))
  ErrorCheck[i, 7]   <- sum(Transecta$N.pts.per.species)
  ErrorCheck[i, 8]   <- sum(Transecta[Transecta$spp.Name == "Unassigned", ]$N.pts.per.species)
  
  ## Agregar cuadratas con multiples appends a la matriz de errores
  if(ErrorCheck$Total_points[i] > 25*ErrorCheck$Quadrats[i]) {
    
    Transect_Split <- split(Transecta, f = Transecta$Image)
    
    ExtraAppend <- NULL
    for(k in 1:length(Transect_Split)) {
      if(sum(Transect_Split[[k]][, 5]) > Transect_Split[[k]][1, 9]) {
        ExtraAppend <- c(ExtraAppend, Transect_Split[[k]][1, 1])
      }
    }
    ErrorCheck[i, 9] <- paste(ExtraAppend, collapse = "; ")
  } else {
    ErrorCheck[i, 9] <- 0
  }
  
  ## Agregar cuadratas que contienen puntos sin asignar
  if(ErrorCheck$Unassigned[i] > 0) {
    ErrorCheck[i, 10]  <- paste(unique(Transecta[Transecta$spp.Name == "Unassigned", 1]), collapse = "; ")
  } else {
    ErrorCheck[i, 10] <- 0
  }
}

ErrorCheck