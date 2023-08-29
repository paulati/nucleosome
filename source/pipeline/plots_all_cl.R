library(ggtext)
library(RColorBrewer)

genome_id <- 'TcruziCL'
source_base_path <- '/home/paula/2023/josefina/nucleosome'
plots_base_path <- '/u01/home/paula/2023/josefina/nucleosome/output/plots/cl'
matrix_base_path <- '/u01/home/paula/2020/nucleosome_github/output/plots'
out_base_path <- '/u03/paula/2021/josefina/MNase'



plots_source_path <- file.path(source_base_path, 'tools/tcruzi_plot2DO/source/plot2DO_plots.R')
file.exists(plots_source_path)
source(plots_source_path)

main <- function() {

  tool <- 'bowtie2' #'hisat2'
  genomes <- c('es', 'es_noEs', 'all')
  params <- list()
  params$colorScaleMax <- 0.05 #valor maximo de la escala, arriba de esto se asigna 4.5
  params$plotType <- "OCC"
  params$reference <- "TSS"
  params$align <- "fivePrime"
  params$simplifyPlot <- FALSE
  params$squeezePlot <- FALSE
  
  params$siteLabel <- "fivePrime"
  
  plots <- list()
  
  for (genome in genomes) {
    
    matrix_file_path <- file.path(matrix_base_path,
                                  paste0('repl1/',
                                         tool, '/clean/', genome,
                                         '/2D_occ_Sites/OCC_matrix.Sites.90_180.CLJ_1_70U_comb.RData'))
    
    file.exists(matrix_file_path)
    
    genes_matrix_file_path <- file.path(matrix_base_path,
                                        paste0('repl1/',
                                               tool, '/clean/', genome,
                                               '/2D_occ_Sites/genes_OCC_matrix.Sites.90_180.CLJ_1_70U_comb.RData'))
    
    file.exists(genes_matrix_file_path)
    
    plots[[genome]][['reads']] <- PlotReadsLength(genome, matrix_file_path, params)
    plots[[genome]][['genes']] <- PlotReadsGenes(genome, genes_matrix_file_path, params)
    
  }
  
  #sacar esto de aca!!!
  graphicalParams <- GetGraphicalParams_custom(params$simplifyPlot, params$squeezePlot)

  plots_in_grid <- list(
                        plots[['es']]$reads$avg_occupancy, 
                        plots[['es']]$reads$heatmap,
                        plots[['es']]$reads$reads_length,
                        plots[['es']]$reads$heatmap_legend,
                        plots[['es']]$genes$heatmap,
                        plots[['es']]$genes$heatmap_legend,
                        plots[['es_noEs']]$reads$avg_occupancy, 
                        plots[['es_noEs']]$reads$heatmap,
                        #plots[['es_noEs']]$reads$reads_length,
                        #plots[['es_noEs']]$reads$heatmap_legend,
                        plots[['es_noEs']]$genes$heatmap,
                        #plots[['es_noEs']]$genes$heatmap_legend,
                        plots[['all']]$reads$avg_occupancy, 
                        plots[['all']]$reads$heatmap,
                        #plots[['all']]$reads$reads_length,
                        #plots[['all']]$reads$heatmap_legend,
                        plots[['all']]$genes$heatmap
                        #plots[['all']]$genes$heatmap_legend
                        )
  
  #https://community.rstudio.com/t/common-axis-title-in-grid-arrange/96353
  # https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html
  # http://www.sthda.com/english/wiki/wiki.php?id_contents=7930
  # https://gist.github.com/tomhopper/faa24797bb44addeba79 <- este
  result.grob <- arrangeGrob(grobs = plots_in_grid,
                                          ncol = ncol(graphicalParams$layout), 
                                          nrow = nrow(graphicalParams$layout),
                                          layout_matrix = graphicalParams$layout,
                                          widths = graphicalParams$gridWidths,
                                          heights = graphicalParams$gridHeights)

  grid.draw(result.grob)
  
  plotFileName <- paste0(tool, "_figure_", genome_id, ".pdf" )
  plotFilePath <- file.path(out_base_path, plotFileName)
  ggsave(plotFilePath, result.grob, width = graphicalParams$plotWidth, height = graphicalParams$plotHeight, 
         units = "in", dpi = 300, scale = 1)   
  
}


PlotReadsLength <- function(genome, matrix_file_path, params) {
  
  # PlotFigure:
  
  load(matrix_file_path)
  
  occMatrix[occMatrix < 0] = 0 # eliminate rounding errors
  if (! is.null(params$colorScaleMax)) {
    occMatrix[occMatrix > params$colorScaleMax] = params$colorScaleMax 
    # override the default colorbar and set everything above the threshold with the threshold value
  }
  
  plotsConfig <- GetPlotsLabels(sampleName, params$plotType, params$reference, params$siteLabel, params$align)
  
  if(genome == "es") {
    plotsConfig$average$mainTitle <- "Esmeraldo"  
  } else if(genome == "es_noEs") {
    plotsConfig$average$mainTitle <- "Es_U_nonEs"  
  } else if(genome == "all") {
    plotsConfig$average$mainTitle <- "CL_Brener_all"  
  }
  
  plotsConfig$fragmentLength$xTitle <-""
  plotsConfig$fragmentLength$mainTitle <- "Es"
  plotsConfig$average$xTitle <- "Position relative to TAS (bp)"
  plotsConfig$average$yTitle <- "Average nucleosome <br/> occupancy"
  plotsConfig$heatmap$xTitle <- "Position relative to TAS (bp)"
  plotsConfig$heatmap$yTitle <- "Fragment length (bp) <br/>  "
  
  heatmapTitles <- plotsConfig$heatmap
  
  graphicalParams <- GetGraphicalParams_custom(params$simplifyPlot, params$squeezePlot)
  
  heatmapTitles$legendTitle <- "Relative coverage (%)"
  
  heatmap <- PlotHeatmap_custom(occMatrix, heatmapTitles$xTitle, heatmapTitles$yTitle, 
                                heatmapTitles$mainTitle, heatmapTitles$legendTitle, 
                                beforeRef, afterRef, lMin, lMax, 
                                graphicalParams$heatmapTheme, graphicalParams$legendTitleTheme,
                                graphicalParams$legendLabelTheme,
                                graphicalParams$scaleXPosition, graphicalParams$scaleYPosition,
                                params$colorScaleMax)
  
  averageTitles <- plotsConfig$average
  avgOccupancy <- PlotAverageOccupancy_custom(occMatrix, beforeRef, afterRef, 
                                       averageTitles$xTitle, averageTitles$yTitle, averageTitles$mainTitle,
                                       graphicalParams$avgOccupancyTheme) 

  fragmentLengthTitles <- plotsConfig$fragmentLength
  fragmentLength <- PlotFragmentLength(lengthHist, lMin, lMax, 
                                       fragmentLengthTitles$xTitle, fragmentLengthTitles$yTitle,
                                       graphicalParams$fragmentLengthTheme)
  
  fragmentLength <- fragmentLength + coord_fixed(ratio = 15) + 
    coord_flip() + theme(aspect.ratio=2/1) 
  
  yBreaks <- c(0.08, 1.5, 3.08)
  yLabels <- c("0.08", "1.5", "3.08")
  yLimits <- c(0, 3.5)
  fragmentLength <- fragmentLength + 
    scale_y_continuous(breaks = yBreaks, labels = yLabels, 
                                      limits = yLimits, expand = c(0,0),
                                          sec.axis = dup_axis())
  
  # Separate legend from heatmap:
  legend_aux <- GetHeatmapLegend(heatmap) 
  heatmapWithoutLegend <- heatmap + theme(legend.position="none")
  heatmapWithoutLegend <- heatmapWithoutLegend + theme(aspect.ratio=1) + coord_fixed() # preserve the aspect ratio
  
  
  result <- list('heatmap' = heatmapWithoutLegend,
                           'heatmap_legend' = legend_aux,
                           'avg_occupancy' = avgOccupancy,
                           'reads_length' = fragmentLength)
                           
  return(result)
  
}
  
PlotReadsGenes <- function(genome, genes_matrix_file_path, params) {

  
  load(genes_matrix_file_path)
  
  params$colorScaleMax <- 3 #valor maximo de la escala, arriba de esto se asigna 4.5
  occMatrix[occMatrix < 0] = 0 # eliminate rounding errors
  
  if (! is.null(params$colorScaleMax)) {
    occMatrix[occMatrix > params$colorScaleMax] = params$colorScaleMax 
    # override the default colorbar and set everything above the threshold with the threshold value
  }
  
  plotsConfig <- GetPlotsLabels(sampleName, params$plotType, params$reference, siteLabel, params$align)
  
  plotsConfig$heatmap$legendTitle <- 'Relative coverage (%)'
  plotsConfig$heatmap$xTitle <- "TAS"
  plotsConfig$heatmap$yTitle <- "Genes <br/>  "
  
  heatmapTitles <- plotsConfig$heatmap
  
  graphicalParams <- GetGraphicalParams_custom(params$simplifyPlot, params$squeezePlot, TRUE)
  
  heatmapTitles$legendTitle  <- 'Relative coverage (%)'
  
  heatmap_genes <- PlotHeatmap_custom(occMatrix, heatmapTitles$xTitle, heatmapTitles$yTitle, 
                                      heatmapTitles$mainTitle, heatmapTitles$legendTitle, 
                                      beforeRef, afterRef, lMin, lMax, 
                                      graphicalParams$heatmapTheme, graphicalParams$legendTitleTheme,
                                      graphicalParams$legendLabelTheme,
                                      graphicalParams$scaleXPosition, graphicalParams$scaleYPosition,
                                      params$colorScaleMax,
                                      TRUE)

  # Separate legend from heatmap:
  legend_aux_genes <- GetHeatmapLegend(heatmap_genes)
  heatmapWithoutLegend_genes <- heatmap_genes + 
    theme(legend.position="none", aspect.ratio=3/2) + coord_fixed() # preserve the aspect ratio
  
  result <- list('heatmap' = heatmapWithoutLegend_genes,
                 'heatmap_legend' = legend_aux_genes,
                 'avg_occupancy' = NULL,
                 'reads_length' = NULL)
  
  return(result)  
    
}

GetHeatmapBreaksAndLabels_custom <- function(occMatrix, colorScaleMax, plotGenes = FALSE) {
  
  if(plotGenes) {
    
    if (is.null(colorScaleMax)) {
      maxValue <- max(occMatrix)
    } else {
      maxValue <- colorScaleMax      
    }
    
    breaks <- round(seq(0, colorScaleMax, length.out = 10), 1)
    labels <- breaks # * 100 # Use percentages instead of fractional numbers  
    limits <- c(0, maxValue)
    
  } else {

    if (is.null(colorScaleMax)) {
      maxValue <- max(occMatrix)
    } else {
      maxValue <- colorScaleMax      
    }
    
    if (maxValue > 0.1) {
      step <- 0.05
    } else if (maxValue > 0.05) {
      step <- 0.01
    } else if (maxValue > 0.01) {
      step <- 0.005
    } else {
      step <- 0.002
    }
    
    breaks <- seq(0, maxValue, step)
    labels <- breaks * 100 # Use percentages instead of fractional numbers  
    limits <- c(0, maxValue)
    
    
  }
  

  result <- list(breaks = breaks, labels = labels, limits = limits)
  
  return(result)
}


PlotHeatmap_custom <- function(occMatrix, xTitle, yTitle, mainTitle, legendTitle,
                        beforeRef, afterRef, lMin, lMax, 
                        customTheme, legendTitleTheme, legendLabelTheme, 
                        scaleXPosition, scaleYPosition, colorScaleMax, plotGenes = FALSE) {
  
  occMatrixMelt <- melt(t(occMatrix))
  
  breaksLabels <- GetHeatmapBreaksAndLabels_custom(occMatrix, colorScaleMax, plotGenes)
  
  # https://rstudio-pubs-static.s3.amazonaws.com/716667_6f5a92af699f4909a91280fdf7ca9bbd.html
  # 100 colors (diverging color scale) with white in the middle:
  getPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
  color_count <- 101
  mycolors <- getPalette(color_count)

  result <- ggplot(occMatrixMelt, aes(x = Var1, y = Var2, fill = value)) + 
    geom_raster(aes(fill = value), interpolate = TRUE) +
    scale_color_gradientn(colors = mycolors, 
                          aesthetics = "fill", 
                          breaks = breaksLabels$breaks, 
                          labels = breaksLabels$labels,
                          limits = breaksLabels$limits
                          )
  
  xBreaks <- seq(-beforeRef, afterRef, by=100) + beforeRef + 1
  if (length(xBreaks) >= 5){
    xLabels <- FixTickLabels(as.character(seq(-beforeRef, afterRef, by=100)), 4, "")  # skip 4 labels
  } else {
    xBreaks <- seq(-beforeRef, afterRef, by=10) + beforeRef + 1
    xLabels <- FixTickLabels(as.character(seq(-beforeRef, afterRef, by=10)), 4, "")  # skip 4 labels
  }
  xLimits <- c(1, beforeRef + afterRef + 1)
  
  if(plotGenes) {

    ymin <- 1
    ymax <- dim(occMatrix)[1]
    
    yBreaks <- seq(ymin, ymax, by=1000) - ymin + 1
    yLabels <- rep("000", length(yBreaks)) 
    yLimits <- c(1, ymax - ymin + 1)
    
  } else {

    yBreaks <- seq(lMin, lMax, by=10) - lMin + 1
    yLabels <- FixTickLabels(as.character(seq(lMin, lMax, by=10)), 4, "")  # skip 4 labels
    yLimits <- c(1, lMax - lMin + 1)

  }
  
  # https://ggplot2.tidyverse.org/reference/guide_colourbar.html
  guideColourbar <- guide_colourbar(title = legendTitle, reverse = FALSE, title.position = "left",
                                    frame.colour = "black", frame.linewidth = 0.5,
                                    ticks.colour = "black", ticks.linewidth = 0.5, 
                                    label.position = "left", 
                                    title.theme = legendTitleTheme,
                                    label.theme = legendLabelTheme,
                                    title.hjust = 0.5,
                                    #barheight = 10,
                                    barheight = 8,                                    
                                    nbin = 1000)  # necessary to shift zero tick to bottom
  
  result <- result + labs(x = xTitle, y = yTitle, title = mainTitle)
  
  scaleX <- scale_x_continuous(breaks = xBreaks, labels = xLabels, limits = xLimits,
                               expand = c(0,0), position = scaleXPosition, 
                               sec.axis = dup_axis())
  
  scaleY <- scale_y_continuous(breaks = yBreaks, labels = yLabels, limits = yLimits,
                               expand = c(0,0), position = scaleYPosition,
                               sec.axis = dup_axis())
  
  result <- result + 
    guides(fill = guideColourbar) +
    scaleY + scaleX + 
    geom_vline(xintercept=beforeRef+1, linetype='longdash', color="white", size=0.4) +
    customTheme
  
  return(result)
  
}  


PlotAverageOccupancy_custom <- function(occMatrix, beforeRef, afterRef, 
                                        xTitle, yTitle, mainTitle, customTheme) {
  
  #setting:
  MAX_AVG_OCC <- 1.2
  
  avgOcc <- colSums(occMatrix)
  avgOcc.df <- as.data.frame(avgOcc)
  avgOcc.df$x <- seq(-beforeRef, afterRef, 1)
  
  deltaYBreaks = max(0.1, round(1.05 * MAX_AVG_OCC/10, digits = 1))
  yBreaks <- seq(0, 1.05 * MAX_AVG_OCC, deltaYBreaks)
  
  if (length(yBreaks) >= 8){
    yLabels <- FixTickLabels(as.character(yBreaks), 1, " ")  # yBreaks
  } else {
    yLabels <- FixTickLabels(as.character(yBreaks), 0, " ")  # yBreaks
  }

  yLimits <- c(0, 1.05 * MAX_AVG_OCC)
  
  
  xBreaks <- seq(-beforeRef, afterRef, 100)
  xLabels <- FixTickLabels(as.character(xBreaks), 4, "") # xBreaks
  xLimits <- c(-beforeRef, afterRef)
  
  xlineInterceps <- seq(-5000, 5000, by=500)
  
  result <- ggplot(data=avgOcc.df, aes(x = x, y = avgOcc)) + 
    geom_line(color="blue") + 
    scale_y_continuous(breaks = yBreaks, labels = yLabels, limits = yLimits, expand = c(0,0),
                       sec.axis = dup_axis()) + 
    scale_x_continuous(breaks = xBreaks, labels = xLabels, limits = xLimits, expand = c(0,0),
                       sec.axis = dup_axis()) + 
    geom_vline(xintercept = xlineInterceps, linetype = 'longdash', color = "lightgray", size = 0.4) +
    labs(x = xTitle, y = yTitle, title = mainTitle) + customTheme
  
  result <- result + coord_fixed(ratio = 1000/1)
  return(result)
  
}



GetGraphicalParams_custom <- function(simplifyPlot, squeezePlot, genesPlot=FALSE) {
  
  baseTheme <- theme_bw() + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5,  vjust = 0.5, colour="black"),
    axis.text.y = element_text(size = 10, angle = 0, hjust = 0.5,  vjust = 0.5, colour="black"),
    axis.title = element_text(size = 12, hjust = 0.5, vjust = 0.5),
    axis.ticks.length = unit(.15, "cm")
  )
  
  legendTitleTheme <- element_text(size = 12, angle = 90, hjust = 0.5, vjust = 0.5, 
                                   margin = margin(t = 0, r = 0, b = 0, l = 0))
  legendLabelTheme <- element_text(size = 10, angle = 0, hjust = 0.5)
  
    layout <- cbind(c(NA, 4, 6), 
                    c(1,2, 5), 
                    c(7, 8, 11),
                    c(13, 14, 17),
                    c(NA,3, NA))  
    gridWidths <- c(1, 1, 1, 1, 1)
    gridHeights <- c(2, 3, 4)
    
    plotHeight <- 9 #12
    plotWidth <- 15 #12
    
    scaleXPosition <- "top"
    scaleYPosition <- "right"
    
    # i have to duplicate axis and add element text to align all the plots.
    # they should use that space but not be shown
    
    if(genesPlot) {
      textColor <- 'white'
    } else {
      textColor <- 'black'
    }
    
    heatmapTheme <- baseTheme + theme(#plot.background = element_rect(fill = "darkblue"),
                                      plot.title = element_blank(),                                                                           
                                      axis.ticks.x.bottom = element_blank(),                                      
                                      axis.ticks.y.left = element_blank(),
                                      axis.text.y.left = element_text(color="white"),                                      
                                      axis.text.y.right = element_markdown(color=textColor, margin = margin(t = 0, r = 0, b = 0, l = 0)),
                                      axis.text.x.bottom = element_text(color="white"),
                                      axis.title.x.bottom = element_text(color="black", margin = margin(t = 5, r = 0, b = 0, l = 0)),
                                      axis.title.x.top = element_blank(),                                      
                                      axis.title.y.left = element_markdown(color="black", margin = margin(t = 0, r = 0, b = 0, l = 0)),
                                      axis.title.y.right = element_blank())
    
    avgOccupancyTheme <- baseTheme + theme(#plot.background = element_rect(fill = "red"),
                                           axis.ticks.x.top = element_blank(),
                                           axis.text.x.top = element_blank(),
                                           axis.title.x.top = element_blank(),
                                           axis.title.x.bottom = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
                                           axis.ticks.y.right = element_blank(),
                                           axis.text.y.right = element_text(color="white"),
                                           axis.title.y.right = element_blank(),
                                           axis.title.y.left = element_markdown(margin = margin(t = 0, r = 0, b = 0, l = 0)))
    
    fragmentLengthTheme <- baseTheme + theme(#plot.background = element_rect(fill = "yellow"),
                                             plot.margin = margin(t = 15, r = 0, b = 15, l = 0),
                                             axis.ticks.x.top = element_blank(),
                                             axis.text.x.top = element_text(color="white"),
                                             axis.title.x.top = element_blank(),                                             
                                             axis.ticks.y.right = element_blank(),
                                             axis.text.y.right = element_text(margin = margin(t = 7, r = 7, b = 7, l = 7)),
                                             axis.title.x.bottom = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
                                             axis.title.y.right = element_blank(),
                                             axis.title.y.left = element_text(margin = margin(t = 7, r = 7, b = 0, l = 7)))
    
  
  result = list(layout = layout, 
                gridWidths = gridWidths, gridHeights = gridHeights, 
                plotWidth = plotWidth, plotHeight = plotHeight, 
                heatmapTheme = heatmapTheme,
                legendTitleTheme = legendTitleTheme,
                legendLabelTheme = legendLabelTheme, 
                avgOccupancyTheme = avgOccupancyTheme,
                fragmentLengthTheme = fragmentLengthTheme,
                scaleXPosition = scaleXPosition,
                scaleYPosition = scaleYPosition)
  
  return(result)
  
}


main()


