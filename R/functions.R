require(Polychrome)

#colors
glasbey <- glasbey.colors(32)
colors_many <- toupper(read_csv(file.path('data','20180313_trubetskoy_colors.csv'))[[3]])
colors <- toupper(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999'))
colors_pat <- toupper(c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5'))
colors_fig <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788",'#984EA3', "light grey", "grey", "dark grey", "#696969")

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

#Marimekko plot without stats
mosaicGG2 <- function(data, X, FILL, colors = colors_many, rect_col = 'white', line_width = 0.25) {
  require(dplyr)
  require(reshape2)
  #require(ggthemes)
  # Proportions in raw data
  DF <- as.data.frame.matrix(table(data[[X]], data[[FILL]]))
  DF$groupSum <- rowSums(DF)
  DF$xmax <- cumsum(DF$groupSum)
  DF$xmin <- DF$xmax - DF$groupSum
  DF$X <- row.names(DF)
  DF$groupSum <- NULL
  DF_melted <- melt(DF, id = c("X", "xmin", "xmax"), variable.name = "FILL")
  DF_melted <- DF_melted %>%
    group_by(X) %>%
    mutate(ymax = cumsum(value/sum(value)),
           ymin = ymax - value/sum(value))
  
  # Chi-sq test
  results <- chisq.test(table(data[[FILL]], data[[X]])) # fill and then x
  resid <- reshape2::melt(results$residuals)
  names(resid) <- c("FILL", "X", "residual")
  
  # Merge data
  DF_all <- merge(DF_melted, resid)
  
  # Positions for labels
  DF_all$xposn <- DF_all$xmin + (DF_all$xmax - DF_all$xmin)/2
  index <- DF_all$xmax == max(DF_all$xmax)
  #DF_all$yposn <- DF_all$ymin[index] + (DF_all$ymax[index] - DF_all$ymin[index])/2
  yposn = 0
  # Plot
  g <- ggplot(DF_all, aes(ymin = ymin,  ymax = ymax, xmin = xmin,
                          xmax = xmax, fill = FILL)) +
    geom_rect(col = rect_col, lwd = line_width) +
    geom_text(aes(x = xposn, label = X),
              y = 1, size = 3, angle = 90, hjust = 1, show.legend = FALSE,check_overlap = T) +
    geom_text(aes(x = max(xmax),  y = yposn, label = FILL),
              size = 3, hjust = 1, show.legend = FALSE,check_overlap = T) +
    scale_fill_manual(FILL, values = colors) +
    scale_x_continuous(X, expand = c(0,0)) +
    scale_y_continuous("Proportion", expand = c(0,0)) +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(g)
}

#hypergeometric enrichment analysis
hyper_test_n <- function(data = df, var1 = "Cluster", var2 = "Region") {
  require(tidyverse) 
  require(broom)
  
  .df <- data.frame()
  for (i in unique(data[[var2]])) {
    data2 <- data
    data2[[var2]] <- factor(ifelse(data2[[var2]] == i, i, paste0("non_",i)), levels = c(i, paste0("non_",i)))
    clusters <- as_tibble(table(data2[[`var1`]]), .name_repair = 'unique')
    colnames(clusters) <- c(var1, 'cluster_size')
    vars <- as_tibble(table(data2[,var1], data2[,var2]), .name_repair = 'unique')
    colnames(vars) <- c(var1, var2, "freq_var2")
    vars_wide <- spread(vars, var2, freq_var2)
    
    vars_df <- vars_wide %>%
      left_join(clusters)
    
    
    #hypergeometric test
    #option a
    test_df<- data.frame(q=vars_df[,i], 
                         m=sum(vars_df[,i]), 
                         n=sum(vars_df[,paste0("non_",i)]),
                         k=vars_df[,4])
    
    colnames(test_df)[1] <- "q"
    
    p_hyper <- apply(test_df, MARGIN = 1, function(x) 1-phyper(max(0,x[[1]]-1), x[[2]], x[[3]], x[[4]])) #probability to get q or more successes in populaton
    
    test_df$p_hyper <- p_hyper
    test_df$Cluster <- vars_df[[`var1`]]
    test_df$enrichment_var <- i
    .df <- .df %>%
      bind_rows(test_df[,c("q","m","n","cluster_size","p_hyper","Cluster","enrichment_var")])
  }
  
  
  .df$padj <- p.adjust(.df$p_hyper, method="BH")
  .df$Significance <- ifelse(.df$padj<0.05 & .df$padj>0.01, '*',
                             ifelse(.df$padj<0.01 & .df$padj>0.001, '**',
                                    ifelse(.df$padj<0.001, '***','n.s.')))
  
  return(.df)
}

test_binom <- function(data = metadata, var1 = "Cluster", var2 = "Cellstate") {
  require(broom)
  .df <- data.frame()
  for (i in unique(data[[var2]])) {
    data2 <- data
    data2[[var2]] <- factor(ifelse(data2[[var2]] == i, i, paste0("non_",i)), levels = c(i, paste0("non_",i)))
    clusters <- as_tibble(table(data2$Cluster), .name_repair = 'unique')
    colnames(clusters) <- c(var1, 'cluster_size')
    vars <- as_tibble(table(data2[,var1], data2[,var2]), .name_repair = 'unique')
    colnames(vars) <- c(var1, var2, "freq_var2")
    vars_wide <- spread(vars, var2, freq_var2)
    
    vars_df <- vars_wide %>%
      left_join(clusters)
    
    
    #binomgeometric test
    #option a
    test_df<- data.frame(q=vars_df[,i], 
                         m=sum(vars_df[,i]), 
                         n=sum(vars_df[,paste0("non_",i)]),
                         k=vars_df[,4])
    colnames(test_df)[1] <- "q"
    p_binom <- apply(test_df, MARGIN = 1, function(x) {glance(binom.test(x = max(0,x[[1]]-1), n = x[[4]], p = x[[2]]/(x[[2]] + x[[3]]), alternative = "greater"))[["p.value"]]}) #probability to get q or more successes in populaton
    
    test_df$p_binom <- p_binom
    test_df$Cluster <- vars_df$Cluster
    test_df$enrichment_var <- i
    .df <- .df %>%
      bind_rows(test_df[,c("q","m","n","cluster_size","p_binom","Cluster","enrichment_var")])
  }
  
  .df$padj <- p.adjust(.df$p_binom, method="BH")
  .df$Significance <- ifelse(.df$padj<0.05 & .df$padj>0.01, '*',
                             ifelse(.df$padj<0.01 & .df$padj>0.001, '**',
                                    ifelse(.df$padj<0.001, '***','n.s.')))
  
  
  
  return(.df)
}

