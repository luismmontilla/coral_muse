randomchooseNpoints <- function(x, nquadrats = 15, npoints = 25) {
  Transecta <- x
  
  if(nquadrats < length(unique(Transecta$Image))) {
    
    Transecta <- Transecta[Transecta$Image %in% sample(unique(Transecta$Image), nquadrats), ]
    
  } else if(nquadrats > length(unique(Transecta$Image))) {
    AddedQuads <- sample(unique(Transecta$Image), nquadrats - length(unique(Transecta$Image)), replace = T)
    
    for(k in 1:length(AddedQuads)) {
      
      NewQuad <- Transecta[Transecta$Image == AddedQuads[k], ]
      
      NewQuad$Image <- paste(NewQuad$Image, "_B", LETTERS[k], sep = "")
      
      Transecta <- rbind(Transecta, NewQuad)
    }
    
  }
  
  Transectasplit <- split(Transecta, f = Transecta$Image)
  
  Transecta100 <- Transectasplit[[1]][0, ]
  Transecta25 <- Transectasplit[[1]][0, ]
  
  for(i in 1:length(Transectasplit)) {
    for(j in 1:length(Transectasplit[[i]]$spp.ID)) {
      
      spoint <- Transectasplit[[i]][j, 5] # Puntos por especie
      
      Transectsplit_i <- Transectasplit[[i]][j, ][rep(seq_len(nrow(Transectasplit[[i]][j, ])), each = spoint),]
      Transectsplit_i$N.pts.per.species <- 1
      Transectsplit_i
      
      Transecta100 <- rbind(Transecta100, Transectsplit_i)
    }
    
    Transecta25 <- rbind(Transecta25, Transecta100[sample(1:100, size = npoints), ]) #En el caso nuevo esta linea cambia y condiciona el numero de puntos
    Transecta100 <- Transectasplit[[1]][0, ]
  }
  
  names(Transecta25) <- names(Transecta)
  
  Transecta <- aggregate(N.pts.per.species ~ ., data = Transecta25, FUN = sum)[, c(1:4, 9, 5:8)]
  Transecta <- Transecta[order(Transecta$Image), ]  
  Transecta
}
