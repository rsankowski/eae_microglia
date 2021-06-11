

#plot expression seurat
plot_expmap_seurat <- function(features, object=all, reduction = "umap", dims=c(1,2), point_size=1, logsc=FALSE, line_width=0, .retain_cl = retain_cl) {
  
  dims <- paste0(Key(object = object[[reduction]]), dims)
  data <- FetchData(object = object, vars = c(dims, "ident", features),  slot = "data")
  
  if (ncol(data) > 4) {
    data2 <- data.frame(data[,1:3], rowSums(data[, 4:ncol(data)]))
  } else {
    data2 <- data
  }
  
  l <- data2[[4]][which(data$ident %in% .retain_cl)]
  mi <- min(l)
  ma <- max(l)
  ColorRamp <- colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)
  ColorLevels <- seq(mi, ma, length = length(ColorRamp))
  v <- round((l - mi)/(ma - mi) * 99 + 1, 0)
  
  kk <- bind_cols(data.frame('l'=l), data[, dims][which(data$ident %in% .retain_cl),]) %>% arrange(l)
  colnames(kk)[2:3] <- c("UMAP_1", "UMAP_2")
  
  if(logsc) {
    plot <- ggplot(kk, aes(UMAP_1, UMAP_2, color = log(l+0.1))) +
      geom_point(size = point_size, pch = 19) +
      scale_color_gradientn('', colors = ColorRamp) +
      theme_void() +
      labs(title = paste(features, collapse = ',')) 
    return(plot)
  }
  else {
    plot <- ggplot(kk, aes(UMAP_1, UMAP_2, color = l)) +
      geom_point(size = point_size, pch = 19) +
      scale_color_gradientn('', colors = ColorRamp) +
      theme_void() +
      labs(title = paste(features, collapse = ','))
    return(plot)
  }
  
}
