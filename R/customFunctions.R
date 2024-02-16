


GenerateContrasts <- function(DESeq_object, 
                              group1 = NULL, 
                              group2 = NULL,
                              meta_data = NULL,
                              mean_norm_counts = NULL,
                              tpc_filtered_list = NULL
                              ){
  
  #print(group1)
  #print(group2)
  
  # if(length(group1) < 1 & length(group2) < 1){
  #   return("Error: No elements were selected to draw contrasts.")
  # }
  # else if (length(group2) >= 1 & length(group1) < 1){
  #   return("Error: Please also select elements in contrast group 1.")
  # }
  # else if (length(group1) >= 1 & length(group2) < 1){
  #   return("Error: Please also select elements in contrast group 2.")
  # }
  # else if(length(intersect(group1, group2)) >1){
  #   return(paste0("Error: The groups to draw the contrasts should not ", 
  #                 "contain any overlaps."))
  # }
  # else if(length(group1) == 1 & length(group2) < 1){
  #   
  #   result <- as.data.frame(results(object = DESeq_object,
  #                                   test = "Wald",
  #                                   name = group1,
  #                                   cooksCutoff = F,
  #                                   independentFiltering = F,
  #                                   parallel = T)
  #                           )
  #   
  #   result$ensembl_gene_id <- gsub(pattern = "\\.[0-9]+$", replacement = "",
  #                            x = rownames(result),
  #                            ignore.case = T)
  #   
  #   result <- merge(x = meta_data, y = result, 
  #                   by = "ensembl_gene_id", all.x = F, all.y = T)
  #   
  #   rownames(result)  <- c(1:nrow(result))
  #   
  #   colnames(result) <- c("Gene.ID", "Gene.Symbol", "Description", 
  #                         "Mean.Counts", "log2.Fold.Change", "Standard.Error", 
  #                         "Statistic", "p.value", "Adjusted.p.value")
  #   
  #   return(result[result$Mean>0,])
  # }
  # # else if (length(group2) < 1 & length(group1) > 1){
  # #   
  # #   result <- as.data.frame(results(object = DESeq_object, 
  # #                                   test = "Wald", 
  # #                                   contrast = list(group1[1], group1[-1]),
  # #                                   listValues = c(1, -1/(length(group1)-1)), 
  # #                                   cooksCutoff = F, 
  # #                                   independentFiltering = F, 
  # #                                   parallel = T))
  # #   
  # #   result$ensembl_gene_id <- gsub(pattern = "\\.[0-9]+$", replacement = "",
  # #                                  x = rownames(result),
  # #                                  ignore.case = T)
  # #   
  # #   result <- merge(x = meta_data, y = result, 
  # #                   by = "ensembl_gene_id", all.x = F, all.y = T)
  # #   
  # #   rownames(result)  <- c(1:nrow(result))
  # #   
  # #   colnames(result) <- c("Gene.ID", "Gene.Symbol", "Description", "Mean", 
  # #                         "log2.Fold.Change", "Standard.Error", "Statistic", 
  # #                         "p.value", "Adjusted.p.value")
  # #   
  # #   return(
  # #     list(paste0("Warning: More than one element was selected in a single ", 
  # #                 "group. The contrast has been drawn between the first ",
  # #                 "element and the remaining elements of the selection."),
  # #          result[result$Mean>0,]
  # #          )
  # #     )
  # # }
  # else{
    
    result <- as.data.frame(results(object = DESeq_object, 
                                    test = "Wald", 
                                    contrast = list(group1, group2),
                                    listValues = c(1/length(group1), 
                                                   -1/length(group2)), 
                                    cooksCutoff = F, 
                                    independentFiltering = F, 
                                    parallel = T))
    
    #calculating actual log2 FC for filtering uncalculable genes
    {
      group1 <- gsub(pattern = "tissue", 
                                   replacement = "", 
                                   x = group1, 
                                   ignore.case = T)
      group1 <- gsub(pattern = "\\.cellTypeNT", 
                                   replacement = "_Tconv", 
                                   x = group1, 
                                   ignore.case = T)
      group1 <- gsub(pattern = "\\.cellTypeT", 
                                   replacement = "_Treg", 
                                   x = group1, 
                                   ignore.case = T)
      group1 <- gsub(pattern = "lymphnode", 
                                   replacement = "LN", 
                                   x =  group1, 
                                   ignore.case = T)
      group1 <- gsub(pattern = "intraepitheliallayer", 
                                   replacement = "IEL", 
                                   x =  group1, 
                                   ignore.case = T)
      group1 <- gsub(pattern = "laminaproprialymphocyte", 
                                   replacement = "LPL", 
                                   x =  group1, 
                                   ignore.case = T)
      group1 <- gsub(pattern = "peyerspatches", 
                                   replacement = "PP", 
                                   x =  group1, 
                                   ignore.case = T)
      
      group2 <- gsub(pattern = "tissue", 
                     replacement = "", 
                     x = group2, 
                     ignore.case = T)
      group2 <- gsub(pattern = "\\.cellTypeNT", 
                     replacement = "_Tconv", 
                     x = group2, 
                     ignore.case = T)
      group2 <- gsub(pattern = "\\.cellTypeT", 
                     replacement = "_Treg", 
                     x = group2, 
                     ignore.case = T)
      group2 <- gsub(pattern = "lymphnode", 
                     replacement = "LN", 
                     x =  group2, 
                     ignore.case = T)
      group2 <- gsub(pattern = "intraepitheliallayer", 
                     replacement = "IEL", 
                     x =  group2, 
                     ignore.case = T)
      group2 <- gsub(pattern = "laminaproprialymphocyte", 
                     replacement = "LPL", 
                     x =  group2, 
                     ignore.case = T)
      group2 <- gsub(pattern = "peyerspatches", 
                     replacement = "PP", 
                     x =  group2, 
                     ignore.case = T)
    }
    #print(group1)
    #print(group2)
    if(length(group1) < 2 ){
      numerator_temp <- as.data.frame(
        mean_norm_counts[, grepl(pattern = group1, 
                                 x = colnames(mean_norm_counts), 
                                 ignore.case = T), drop = F])
      colnames(numerator_temp) <- group1
    }
    else {
      numerator_temp <- as.data.frame(
        mean_norm_counts[, grepl(pattern = group1[1], 
                                 x = colnames(mean_norm_counts), 
                                 ignore.case = T), drop = F])
      colnames(numerator_temp) <- group1[1]
      for(i in seq(2, length(group1))){
        temp <- as.data.frame(
          mean_norm_counts[, grepl(pattern = group1[i], 
                                   x = colnames(mean_norm_counts), 
                                   ignore.case = T), drop = F])
        numerator_temp <- cbind(numerator_temp, temp)
        colnames(numerator_temp)[i] <- group1[i]
        rm(temp)
      }
    }
    rm(i)
    if(length(group2) < 2 ){
      denominator_temp <- as.data.frame(
        mean_norm_counts[, grepl(pattern = group2, 
                                 x = colnames(mean_norm_counts), 
                                 ignore.case = T), drop = F])
      colnames(denominator_temp) <- group2
    }
    else {
      denominator_temp <- as.data.frame(
        mean_norm_counts[, grepl(pattern = group2[1], 
                                 x = colnames(mean_norm_counts), 
                                 ignore.case = T), drop = F])
      colnames(denominator_temp) <- group2[1]
      for(i in seq(2, length(group2))){
        temp <- as.data.frame(
          mean_norm_counts[, grepl(pattern = group2[i], 
                                   x = colnames(mean_norm_counts), 
                                   ignore.case = T), drop = F])
        denominator_temp <- cbind(denominator_temp, temp)
        colnames(denominator_temp)[i] <- group2[i]
        rm(temp)
      }
    }
    rm(i)
    result$Actuallog2FC <- rowMeans(numerator_temp, na.rm = T) - 
      rowMeans(denominator_temp, na.rm = T)
    result$Actuallog2FC[!is.finite(result$Actuallog2FC)] <- NA
    
    #exclude genes with NA/NaN/Inf/-Inf 
    #result <- result[complete.cases(result),]
    #print(unique(result[!complete.cases(result),]$Actuallog2FC)) ## all NaN
    #print(unique(result[!is.finite(result$Actuallog2FC),]$Actuallog2FC)) ## still all NaN
    #result$Actuallog2FC <- NULL
    
    result$ensembl_gene_id <- gsub(pattern = "\\.[0-9]+$", replacement = "",
                                   x = rownames(result),
                                   ignore.case = T)
    
    result <- merge(x = meta_data, y = result, 
                    by = "ensembl_gene_id", all.x = F, all.y = T)
    
    rownames(result)  <- c(1:nrow(result))
    result <- result[,c(1,2,3,4,10,5,8,9,6,7)]
    # colnames(result) <- c("Gene.ID", "Gene.Symbol", "Description",
    #                       "Mean.Counts", "log2.Fold.Change", "Standard.Error",
    #                       "Statistic", "p.value", "Adjusted.p.value")#,
    #                       #"Actual.log2.Fold.Change")
    colnames(result) <- c("Gene.ID", "Gene.Symbol", "Description",
                          "Mean.Counts", "log2.Fold.Change",
                          "DESeq2.Beta.Coefficients", "p.value",
                          "Adjusted.p.value", "Standard.Error", "Statistic")#,
    #"Actual.log2.Fold.Change")
    
    if(!is.null(tpc_filtered_list)){
      result <- result[result$Gene.ID %in% tpc_filtered_list,]
    }
    
    #return(result[result$Mean>0,])
    #return(result[result$Mean.Counts>0,])
    return(result)
    
  # }
  
  
}


GenerateVolcanoPlot <- function(plot_data,
                                names_group1,
                                names_group2,
                                highlight = NULL, 
                                format = "svg", 
                                width = 880, 
                                height = 562,
                                log2FCFilt = NULL,
                                AdjpFilt = NULL
                                ){
  
  #checking if the results were accompanied by a warning or error text
  if(class(plot_data) == "list"){
    plot_data <- plot_data[[2]]
  }
  
  #assigning default values in case they were not forwarded erroneously
  if(is.null(log2FCFilt)) log2FCFilt <- 2
  if(is.null(AdjpFilt)) AdjpFilt <- 0.01
  
  plot_data$Row.Order <- c(1:nrow(plot_data))
  
  #plot_data <- plot_data[order(plot_data$log2.Fold.Change, decreasing = T),]
  plot_data <- plot_data[order(plot_data$DESeq2.Beta.Coefficients, decreasing = T),]
  #print(summary(plot_data))
  
  #assigning Facet groups to genes with missing/non-calculable Actuallog2FC
  plot_data$Facet <- "Both"
  plot_data$Facet[!is.finite(plot_data$log2.Fold.Change) & 
                    plot_data$DESeq2.Beta.Coefficients > 0] <- "Right"
  plot_data$Facet[!is.finite(plot_data$log2.Fold.Change) & 
                    plot_data$DESeq2.Beta.Coefficients < 0] <- "Left"
  plot_data$Facet <- factor(x = plot_data$Facet, 
                            levels = c("Left", "Both", "Right"))
  
  #assigning the 'max+' log2 fold changes of "Left" and "Right" genes
  nodetect_left <- -(ceiling(max(abs(
    plot_data$DESeq2.Beta.Coefficients[plot_data$Facet == "Both"]))) + 5) #* -1
  nodetect_right <- (ceiling(max(abs(
    plot_data$DESeq2.Beta.Coefficients[plot_data$Facet == "Both"]))) + 5)
  
  plot_data$DESeq2.Beta.Coefficients[plot_data$Facet == "Left"] <- 
    nodetect_left
  plot_data$DESeq2.Beta.Coefficients[plot_data$Facet == "Right"] <- 
    nodetect_right
  
  #formatting plot title
  # if(length(names_group2) < 1){
  #   title <- paste0(paste(names_group1, collapse = ", "),
  #                   "  vs  ", "rest")
  # }
  # else{
  #   title <- paste0(paste(names_group1, collapse = ", "),
  #                   "  vs  ",
  #                   paste(names_group2, collapse = ", "))
  # }
  title <- paste0(paste(names_group1, collapse = ", "), 
                  "  vs  ", 
                  paste(names_group2, collapse = ", "))
  #formatting plot title here instead of in the already bulky plot code
  {
    #title <- names(merged_results_list[[i]][[j]])[k]
    title <- gsub(pattern = "blood", 
                  replacement = "Blood", 
                  x = title, ignore.case = T)
    title <- gsub(pattern = "spleen", 
                  replacement = "Spleen", 
                  x = title, ignore.case = T)
    title <- gsub(pattern = "lymphnode", 
                  replacement = "LN", 
                  x = title, ignore.case = T)
    title <- gsub(pattern = "peyerspatches", 
                  replacement = "PP", 
                  x = title, ignore.case = F)
    title <- gsub(pattern = "intraepitheliallayer", 
                  replacement = "IEL", 
                  x = title, ignore.case = F)
    title <- gsub(pattern = "laminaproprialymphocyte", 
                  replacement = "LPL", 
                  x = title, ignore.case = F)
    title <- gsub(pattern = "kidney", 
                  replacement = "Kidney", 
                  x = title, ignore.case = F)
    title <- gsub(pattern = "lung", 
                  replacement = "Lung", 
                  x = title, ignore.case = F)
    title <- gsub(pattern = "liver", 
                  replacement = "Liver", 
                  x = title, ignore.case = F)
    title <- gsub(pattern = "pancreas", 
                  replacement = "Pancreas", 
                  x = title, ignore.case = F)
    
    # title <- paste0("Comparing ", gsub(pattern = "_", 
    #                                    replacement = ", ", 
    #                                    x = gsub(pattern = "__Vs__", 
    #                                             replacement = "  vs  ", 
    #                                             x = title, 
    #                                             ignore.case = T), 
    #                                    ignore.case = T))
  }
  
  #assign colours based on significance
  plot_data$coloured <- "not_selected"
  # plot_data$coloured[which(plot_data$Adjusted.p.value <= 0.01 &
  #                            (plot_data$log2.Fold.Change >= 2 |
  #                               plot_data$log2.Fold.Change <= -2)
  #                          )] <- "selected"
  plot_data$coloured[which(plot_data$Adjusted.p.value <= AdjpFilt &
                             (plot_data$DESeq2.Beta.Coefficients >= log2FCFilt |
                                plot_data$DESeq2.Beta.Coefficients <= -log2FCFilt)
                           )] <- "selected"
  # plot_data$coloured[grepl(pattern = genes_of_int, 
  #                          x = plot_data$MGI.symbol, 
  #                          ignore.case = T)] <- "selected"
  plot_data <- plot_data[order(plot_data$coloured),]
  
  #print(head(plot_data))
  
  #generating the plot
  scatter_plot <- ggplot(data = plot_data, 
                         aes_(#x = plot_data$log2.Fold.Change, 
                           x = plot_data$DESeq2.Beta.Coefficients, 
                              y = -log10(plot_data$Adjusted.p.value), 
                              #Gene = str_extract(plot_data$Gene.Symbol, '\\w*'),
                              Gene = gsub(pattern = ";.*", 
                                          replacement = "", 
                                          x = plot_data$Gene.Symbol, 
                                          ignore.case = T),
                              #log2_FC = plot_data$log2.Fold.Change,
                           log2_FC = plot_data$log2.Fold.Change,
                           Beta = plot_data$DESeq2.Beta.Coefficients, 
                              Sig = plot_data$Adjusted.p.value,
                              colour = plot_data$coloured, 
                              #key = plot_data$Gene.ID
                              key = plot_data$Row.Order
                              ),
                         ) + 
    # geom_hline(yintercept = 0, colour = "#4d4d4d", size = 0.25) + 
    # geom_hline(yintercept = -log10(0.01), colour = "#4d4d4d", size = 0.25) + 
    # geom_vline(xintercept = 0, colour = "#4d4d4d", size = 0.25) + 
    # geom_vline(xintercept = 2, colour = "#4d4d4d", size = 0.25) + 
    # geom_vline(xintercept = -2, colour = "#4d4d4d", size = 0.25) + 
    geom_point(shape = 19, size = 0.8) + 
    scale_color_manual(labels = c("Other genes", "Genes of\ninterest"),
                       #values = c("grey", "darkviolet")) + 
                       #values = c("grey", "#FF9999")) + 
                       #values = c("grey", "#dd4b24")) + 
                       values = c("darkgrey", "#fd7239")) + 
    labs(colour = "", 
         y = "-log10 Adjusted p.value", 
         x = "DESeq2 Beta Coefficients (log2 scale)", 
         #title = title
         caption = title
         ) + 
    # geom_hline(yintercept = -log10(0.01), linetype="dashed", 
    #            color = "darkgrey") + 
    # geom_vline(xintercept = -5, linetype="dashed", color = "darkgrey") + 
    # geom_vline(xintercept = 5, linetype="dashed", color = "darkgrey") + 
    #scale_fill_discrete(name = "Colour", labels = c("", "1% by quartile")) + 
    #scale_colour_viridis_c(option = "plasma") + 
    theme_light() + 
    #theme_light(base_size = 12) #, legend.position='none') + 
    theme(plot.caption = element_text(hjust = 0, face= "italic"), 
          legend.key.width = unit(0.15,"line")
          )
  
  #temp <- ggplot_build(scatter_plot)
  #print(temp)
  
  # xbreaks <- sort(c(
  #   na.omit(ggplot_build(scatter_plot)$layout$panel_params[[1]]$x$breaks),
  #   -5, 5))
  #xbreaks <- sort(c(0,-2, 2))
  xbreaks <- sort(c(0,-log2FCFilt, log2FCFilt, 5, -5, 10, -10, 15, -15, 20, -20))
  # xbreaks <- 
  #   sort(c(-log2FCFilt, log2FCFilt, 
  #          seq(0, -max(abs(
  #            plot_data$DESeq2.Beta.Coefficients[plot_data$Facet == "Both"])), -5),
  #          min(plot_data$DESeq2.Beta.Coefficients), 
  #          seq(0, max(abs(
  #            plot_data$DESeq2.Beta.Coefficients[plot_data$Facet == "Both"])), 5), 
  #          max(plot_data$DESeq2.Beta.Coefficients)
  #          ))
  xbreaks <- sort(c(-log2FCFilt, log2FCFilt, 
                    seq(0, nodetect_left, -5), nodetect_left, 
                    seq(0, nodetect_right, 5), nodetect_right, 
                    max(plot_data$DESeq2.Beta.Coefficients)
                    )
                  )
  xbreaks <- sort(unique(xbreaks))
  
  ybreaks <- sort(c(
    na.omit(ggplot_build(scatter_plot)$layout$panel_params[[1]]$y$breaks),
    -log10(AdjpFilt)))
  #ybreaks <- sort(c(0, -log10(0.01)))
  
  xlabels <- as.character(xbreaks)
  xlabels[xlabels == min(plot_data$DESeq2.Beta.Coefficients) | 
            xlabels == max(plot_data$DESeq2.Beta.Coefficients)] <- 
    sprintf("Detected\nonly in\nthis group")

  ylabels <- as.character(ybreaks)
  #ylabels[ylabels == -log10(0.01)] <- "\nadjusted\np.value  \n0.01      \n\n\n"
  ylabels[ylabels == -log10(AdjpFilt)] <- sprintf("\nadj.p\n%s\n\n\n", AdjpFilt)

  scatter_plot <- scatter_plot +
    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    scale_y_continuous(breaks = ybreaks, labels = ylabels) + 
    # expand_limits(x = c(-ceiling(max(abs(plot_data$log2.Fold.Change))),
    #                     ceiling(max(abs(plot_data$log2.Fold.Change))))
    #               ) +
    expand_limits(x = c(-ceiling(max(abs(plot_data$DESeq2.Beta.Coefficients))),
                        ceiling(max(abs(plot_data$DESeq2.Beta.Coefficients))))
                  ) +
    theme(panel.grid.minor.x = element_blank(), 
          #panel.grid.major.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          #panel.grid.major.y = element_blank(),
          #axis.line = element_line(colour = 'black')
          panel.grid.major = element_line(colour = "#4d4d4d"),
          axis.ticks = element_line(colour = "#4d4d4d")
          ) + 
    ggExtra::removeGrid()
  
  #adding floating text annotations on the plot area
  #print(title)
  label_LHS <- gsub(pattern = ", ", 
                    replacement = "\n", 
                    x = gsub(pattern = "Comparing ", 
                             replacement = "", 
                             x = unlist(strsplit(x = title, 
                                                 split = "  vs  "))[2], 
                             ignore.case = T), 
                    ignore.case = T)
  label_LHS <- paste0("Expressed higher in,\n", label_LHS)
  label_LHS <- gsub(pattern = "_", 
                    replacement = " ", 
                    x = label_LHS, 
                    ignore.case = T)
  #print(label_LHS)
  label_RHS <- gsub(pattern = ", ", 
                    replacement = "\n", 
                    x = gsub(pattern = "Comparing ", 
                             replacement = "", 
                             x = unlist(strsplit(x = title, 
                                                 split = "  vs  "))[1], 
                             ignore.case = T), 
                    ignore.case = T)
  label_RHS <- paste0("Expressed higher in,\n", label_RHS)
  label_RHS <- gsub(pattern = "_", 
                    replacement = " ", 
                    x = label_RHS, 
                    ignore.case = T)
  
  # scatter_plot <- scatter_plot +
  #   annotate(geom = "text",
  #            x = -ceiling(max(abs(plot_data$log2.Fold.Change))),
  #            y = ceiling(max(-log10(plot_data$Adjusted.p.value))),
  #            label = label_LHS,
  #            size = 2,
  #            hjust = 0,
  #            vjust = 1) +
  #   annotate(geom = "text",
  #            x = ceiling(max(abs(plot_data$log2.Fold.Change))),
  #            y = ceiling(max(-log10(plot_data$Adjusted.p.value))),
  #            label = label_RHS,
  #            size = 2,
  #            hjust = 1,
  #            vjust = 1)
  
  #print(ggplot_build(scatter_plot)$layout$panel_scales_x[[1]]$range$range)
  #print(ggplot_build(scatter_plot)$layout$panel_params)
  ##print(ggplot_build(scatter_plot)$layout$panel_params[[1]]$x$breaks)
  ##print(ggplot_build(scatter_plot)$layout$panel_params[[1]]$y$breaks)
  #print(ggplot_build(scatter_plot)$layout$panel_scales_y[[1]]$range$range)
  #print(ggplot_build(scatter_plot)$layout$panel_params[[1]]$x.major_source)
  
  plot_object <- ggplotly(scatter_plot, 
                  tooltip = c("Gene", "log2_FC", "Sig"), 
                  source = "volplot" 
                  ) %>% 
           hide_legend() %>% toWebGL() %>% 
           onRender(js, data = "TraceMapping") %>% 
           config(displaylogo = F, 
                  displayModeBar = T, 
                  toImageButtonOptions = list(format = format, 
                                              width = width, 
                                              height = height),  
                  modeBarButtonsToRemove = list("sendDataToCloud", "zoom2d", 
                                                "zoomIn2d", "zoomOut2d", 
                                                "pan2d", "select2d", "lasso2d", 
                                                "autoScale2d", 
                                                "hoverClosestCartesian", 
                                                "hoverCompareCartesian"
                                                )
                  ) %>% 
           layout(#title = title, 
                  #titlefont = list(size = 10), 
                  margin = list(l = 10#,
                                #r = 10,
                                #b = 20#,
                                #t = 10,
                                #pad = 4
                                ),
                  annotations = list(
                    list(#x = -ceiling(max(abs(plot_data$log2.Fold.Change))), 
                      x = -ceiling(max(abs(plot_data$DESeq2.Beta.Coefficients))), 
                         y = ceiling(max(-log10(plot_data$Adjusted.p.value))),
                         align = "left",
                         #xref = "paper", yref = "paper",
                         xanchor = "left", 
                         yanchor = "top",
                         #xshift = 0, yshift = 0,
                         text = label_LHS, 
                         font = list(size = 11.5),
                         showarrow = FALSE
                         ), 
                    list(#x = ceiling(max(abs(plot_data$log2.Fold.Change))), 
                      x = ceiling(max(abs(plot_data$DESeq2.Beta.Coefficients))), 
                         y = ceiling(max(-log10(plot_data$Adjusted.p.value))),
                         align = "right",
                         #xref = "paper", yref = "paper",
                         xanchor = "right", 
                         yanchor = "top",
                         #xshift = 0, yshift = 0,
                         text = label_RHS, 
                         font = list(size = 11.5),
                         showarrow = FALSE
                         )
                    )
                  )
  plot_object$x$layout$shapes <- append(plot_object$x$layout$shapes,
                                        list(list(type = "line",
                                                  layer = "below",
                                                  x0 = 0, x1 = 1,
                                                  xref = "paper",
                                                  y0 = 0, y1 = 0,
                                                  line = list(color = "#4d4d4d",
                                                              width = 1)),
                                             list(type = "line",
                                                  layer = "below",
                                                  x0 = 0, x1 = 1,
                                                  xref = "paper",
                                                  y0 = -log10(AdjpFilt), 
                                                  y1 = -log10(AdjpFilt),
                                                  line = list(color = "#4d4d4d",
                                                              width = 1)),
                                             list(type = "line",
                                                  layer = "below",
                                                  y0 = 0, y1 = 1,
                                                  yref = "paper",
                                                  x0 = 0, x1 = 0,
                                                  line = list(color = "#4d4d4d",
                                                              width = 1)),
                                             list(type = "line",
                                                  layer = "below",
                                                  y0 = 0, y1 = 1,
                                                  yref = "paper",
                                                  x0 = log2FCFilt, x1 = log2FCFilt,
                                                  line = list(color = "#4d4d4d",
                                                              width = 1)),
                                             list(type = "line",
                                                  layer = "below",
                                                  y0 = 0, y1 = 1,
                                                  yref = "paper",
                                                  x0 = -log2FCFilt, x1 = -log2FCFilt,
                                                  line = list(color = "#4d4d4d",
                                                              width = 1))
                                             )
                                        )
  
  return(plot_object)
}


GenerateColPlot <- function(mean_norm_counts,
                            selection, 
                            format = "svg", 
                            width = 880, 
                            height = 562
                            ){
  
  #print(selection)
  #selecting the data
  plot_data <- mean_norm_counts[mean_norm_counts$Gene.stable.ID %in% 
                                  selection$Gene.ID,]
  
  
  #formatting plot_data to long format for plotting
  #rownames(plot_data) <- plot_data[,1]
  rownames(plot_data) <- plot_data$Gene.stable.ID
  #plot_data <- plot_data[,-1]
  plot_data$Gene.stable.ID <- NULL
  plot_data <- t(plot_data)
  
  
  temp <- data.frame(Exprs = plot_data[,1], GeneID = colnames(plot_data)[1], 
                     Sample = rownames(plot_data))
  
  if(ncol(plot_data)>1){
    for(i in seq(2, ncol(plot_data))){
      temp <- rbind(temp, data.frame(Exprs = plot_data[,i], 
                                     GeneID = colnames(plot_data)[i], 
                                     Sample = rownames(plot_data), 
                                     stringsAsFactors = F))
    }
  }
  
  plot_data <- as.data.frame(temp)
  
  plot_data <- merge(x = plot_data, 
                     y = selection, 
                     by.x = "GeneID", 
                     by.y = "Gene.ID", 
                     all.x = T, 
                     all.y = F)
  
  #print(head(plot_data))
  
  colnames(plot_data) <- gsub(pattern = "Gene.Symbol", 
                              replacement = "Gene", 
                              x = colnames(plot_data), 
                              ignore.case = T)
  
  plot_data$GeneID <- as.factor(plot_data$GeneID)
  plot_data$Sample <- as.factor(plot_data$Sample)
  #plot_data$Gene <- as.character(str_extract(plot_data$Gene, '\\w*'))
  plot_data$Gene <- as.character(gsub(pattern = ";.*", 
                                      replacement = "", 
                                      x = plot_data$Gene, 
                                      ignore.case = T))
  
  #if gene symbol is missing, substitute with ENMUSG
  plot_data$Gene[which(plot_data$Gene == "" | is.na(plot_data$Gene))] <-
    as.character(plot_data$GeneID[which(plot_data$Gene == "" | 
                                          is.na(plot_data$Gene))])
  #print(plot_data$Gene[which(plot_data$Gene == "")])
  #print(plot_data$GeneID[which(plot_data$Gene == "")])
  
  plot_data$Gene <- as.factor(plot_data$Gene)
  
  plot_data$Sample <- as.factor(gsub(pattern = "_", 
                                     replacement = " ", 
                                     x = plot_data$Sample, 
                                     ignore.case = T))

  #reordering the samples to group by cell type first
  plot_data$Sample <- factor(x = plot_data$Sample, 
                             levels = levels(plot_data$Sample)[
                               match(unique(c(
                                 as.character(
                                   plot_data$Sample[grepl(pattern = "Tconv$",
                                                          x = plot_data$Sample,
                                                          ignore.case = T)]),
                                 as.character(
                                   plot_data$Sample[grepl(pattern = "Treg$",
                                                          x = plot_data$Sample,
                                                          ignore.case = T)])
                                 )), 
                                 plot_data$Sample
                                 )]
                               )
  
  col_plot <- ggplot(data = plot_data, 
                     aes(x = Sample, 
                         #y = log2(Count), 
                         #y = log10(Count), 
                         y = Exprs, 
                         name = GeneID)#,
                         #Counts = Count)
                         #Exprs = Exprs)
                     ) + 
    geom_point(size = 1, aes(colour = Gene)) + #, shape = Source)) + 
    scale_color_viridis(discrete = T) + 
    labs(x = "Tissue and cell type", 
         y = "Mean log2 normalised expression") + 
    #ylim(0, ceiling(max(plot_data$Exprs))) + 
    ylim(ifelse(min(plot_data$Exprs, na.rm = T) < 0, 
                floor(min(plot_data$Exprs, na.rm = T)), 
                0), 
         ceiling(max(plot_data$Exprs))) + 
    theme_light() + 
    theme(legend.title = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1))
    #theme_dark() + 
    # theme(legend.position = "bottom", 
    #       )
  
  return(ggplotly(col_plot, 
                  #tooltip = c("Gene", "Sample", "Counts"), 
                  tooltip = c("Gene", "Sample", "Exprs"), 
                  source = "colplot") %>% toWebGL() %>% 
           config(displaylogo = F, 
                  displayModeBar = T, 
                  toImageButtonOptions = list(format = format, 
                                              width = width, 
                                              height = height), 
                  modeBarButtonsToRemove = list("sendDataToCloud", "zoom2d", 
                                                "zoomIn2d", "zoomOut2d", 
                                                "pan2d", "select2d", "lasso2d", 
                                                "autoScale2d", 
                                                #"hoverClosestCartesian", 
                                                #"hoverCompareCartesian", 
                                                "toggleSpikelines"
                  )
           ) %>% 
           add_annotations(text = "Gene", xref = "paper", yref = "paper", 
                           x = 1.04, 
                           xanchor = "left", 
                           y = 0.9, 
                           yanchor = "bottom",    # Same y as legend below
                           legendtitle = TRUE, showarrow = FALSE ) %>%
           layout(#legend = list(orientation = "h", y = -2), 
                  legend = list(y = 0.9,
                                yanchor="top")#, 
                  #xaxis = list(fixedrange = T), 
                  #yaxis = list(fixedrange = F)
                  )
         )
  
}


TextWrap <- function(name_list,
                     max_char_per_line = 70
                     ){
  
  formatted_text <- ""
  current_text <- ""
  
  for(i in seq(1, length(name_list))){
    
    if(sum(nchar(current_text), nchar(", "), nchar(name_list[i])) > 
       max_char_per_line){
      formatted_text <- paste0(formatted_text, current_text, "\n")
      current_text <- name_list[i]
    }
    else{
      current_text <- paste0(current_text, ", ", name_list[i])
    }
    
    #print(current_text)
    
  }
  
  formatted_text <- paste0(formatted_text, current_text)
  
  return(gsub(pattern = "^, *", 
              replacement = "", 
              x = formatted_text, 
              ignore.case = T))
  
}

TextWrapHTML <- function(name_list,
                     max_char_per_line = 70
                     ){
  
  formatted_text <- ""
  current_text <- ""
  
  for(i in seq(1, length(name_list))){
    
    if(sum(nchar(current_text), nchar(", "), nchar(name_list[i])) > 
       max_char_per_line){
      formatted_text <- paste0(formatted_text, current_text, "<br>")
      current_text <- name_list[i]
    }
    else{
      current_text <- paste0(current_text, ", ", name_list[i])
    }
    
    #print(current_text)
    
  }
  
  formatted_text <- paste0(formatted_text, current_text)
  
  return(gsub(pattern = "^, *", 
              replacement = "", 
              x = formatted_text, 
              ignore.case = T))
  
}


#code for creating deseq object
{
  # 
  # #importing and formatting counts matrix
  # {
  #   counts_matrix <- read.table(file = paste(loc_counts, "counts.txt", sep = ""),
  #                               sep = "\t",
  #                               stringsAsFactors = F,
  #                               header = T)
  #   rownames(counts_matrix) <- counts_matrix[,1]
  #   counts_matrix <- counts_matrix[,-1]
  #   
  #   meta_data <- counts_matrix [,1:5]
  #   counts_matrix <- counts_matrix[,-(1:5)]
  #   
  #   ids <- c(1:ncol(counts_matrix))
  #   counts_matrix[,ids] <- as.numeric( unlist( counts_matrix[,ids] ) )
  #   
  #   # #removing contaminated sample Pancreas_T_3
  #   # counts_matrix <- counts_matrix[,!grepl(pattern = "Pancreas_T_3", 
  #   #                                        x = colnames(counts_matrix), 
  #   #                                        ignore.case = T)]
  #   
  #   counts_matrix <- as.matrix(counts_matrix)
  #   
  # }
  # 
  # #reading in DESeq2 outputs
  # {
  #   #reading norm counts
  #   norm_counts <- read.delim(file = paste0(loc_counts, "normalised_counts.txt"),
  #                             sep = "\t",
  #                             stringsAsFactors = F)
  #   #norm_counts <- norm_counts[rowSums(as.matrix(norm_counts))>0,]
  #   
  #   #reading the DESeq2 analysis output table
  #   merged_results_list <- readRDS(file = paste0(loc_diffexprs, 
  #                                                "merged_results_list.rds"))
  #   
  #   contrast_groups <- list(c("Blood"), 
  #                           c("Spleen", "Lymph Node", "Peyer's Patches"), 
  #                           c("Intraepithelial Layer", 
  #                             "Lamina Propria Lymphocyte"), 
  #                           c("Kidney", "Lung", "Liver", "Pancreas"),
  #                           c("cellTypeT"),
  #                           c("cellTypeNT"))
  #   
  #   #ensuring that names in contrast_groups list are all lower letter and do not 
  #   # contain spaces or punctuation
  #   for(i in seq(1, length(contrast_groups)))
  #   {
  #     contrast_groups[[i]] <- tolower(unique(
  #       gsub("_.*|[[:space:]]|[[:punct:]]", "", contrast_groups[[i]])
  #     ))
  #   }
  #   rm(i)
  #   
  #   output_table <- read.delim(file = paste0(loc_diffexprs_merged, 
  #                                            "output_table.txt"), 
  #                              sep = "\t", stringsAsFactors = F)
  #   
  #   colnames(output_table) <- output_table[3,]
  #   output_table <- output_table[-(1:3),]
  #   
  # }
  # 
  # #additional genes of interest that need to be highlighted
  # genes_of_int <- c("CD44", "CD62L", "CD69", "CD103", "ICOS", "KLRG1", "NK1.1", 
  #                   "NRP1", "ST2", "PD1", "NK1")
  # genes_of_int <- paste0("\\b(", paste(genes_of_int, collapse = "|"), ")\\b")
  # genes_of_int <- gsub(pattern = "\\.", replacement = "\\\\.", x = genes_of_int)
  # 
  # #incomplete code that would be used for the webtool
  # {
  #   #initialising DESeq object
  #   print("initialising DESeq object")
  #   DESeq_object <- DESeqDataSetFromMatrix(countData = counts_matrix, 
  #                                          colData = annotations, 
  #                                          design = ~0+tissue:cellType)
  #   
  #   #adding meta_data to the dataset
  #   print("adding meta_data to the dataset")
  #   mcols(DESeq_object) <- DataFrame(mcols(DESeq_object),  meta_data)
  #   
  #   #prefiltering low/no count genes
  #   #print("prefiltering low/no count genes")
  #   #keep <- rowSums(counts(DESeq_object)) >=10
  #   #DESeq_object <- DESeq_object[keep,]
  #   
  #   #performing the DESeq analysis
  #   print("Performing the default DESeq analysis pipeline")
  #   DESeq_object <- DESeq(object = DESeq_object,
  #                         test = "Wald",
  #                         betaPrior = F, 
  #                         minReplicatesForReplace = Inf, 
  #                         parallel = T)
  #   
  #   result_names <- resultsNames(DESeq_object)
  #   
  #   saveRDS(object = DESeq_object, file = "DESeq_Object.rds")
  #   
  #   # #in case we need to automate the name generation for different formulae
  #   # gsub(pattern = "\\W", replacement = "_", x = design(DESeq_object), ignore.case = T)
  #   # 
  #   # str_split(string = gsub(pattern = "\\W", 
  #   #                         replacement = "_", 
  #   #                         x = design(DESeq_object), 
  #   #                         ignore.case = T), 
  #   #           pattern = "_")
  #   
  #   
  #   names(result_names) <- gsub(pattern = "^_", 
  #                               replacement = "", 
  #                               x = gsub(pattern = "(tissue|cellType|\\.)+", 
  #                                        replacement = "_", 
  #                                        x = result_names, 
  #                                        ignore.case = T), 
  #                               ignore.case = T)
  #   
  #   
  #   #creating a list to store results
  #   results_list <- list()
  #   
  # }
  
}

