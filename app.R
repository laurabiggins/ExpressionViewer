

# Developed by Samar H.K. Tareen


#importing libraries & functions, reading in data, and setting up global options
{
  #setting global options
  options(scipen = 5, digits = 4)
  
  #libraries
  {
    library(ggExtra)
    library(shiny)
    library(DESeq2)
    library(shinyjs)
    #library(shinyBS)
    library(DT)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(viridis)
    library(plotly)
    library(stringr)
    library(htmlwidgets)
    library(shinyWidgets)
    library(shinycssloaders)
    library(png)
    library(htmltools)
    #library(shinythemes)
    #library(shinydashboard)
    #library(dashboardthemes)
    library(rintrojs)
    library(shinyalert)
    
  }
  
  #data
  {
    #deseq2 object
    DESeq_object <- readRDS(file = "data/DESeq_Object.rds")
    
    #formatting names for available contrasts
    results_names <- resultsNames(DESeq_object)
    #print(results_names)
    
    #aesthetic formatting of name labels
    {
      names(results_names) <- gsub(pattern = "tissue", 
                                   replacement = "", 
                                   x = results_names, 
                                   ignore.case = T)
      
      names(results_names) <- gsub(pattern = "\\.cellTypeNT", 
                                   replacement = "_Tconv", 
                                   x = names(results_names), 
                                   ignore.case = T)
      
      names(results_names) <- gsub(pattern = "\\.cellTypeT", 
                                   replacement = "_Tregs", 
                                   x = names(results_names), 
                                   ignore.case = T)
      
      # names(results_names) <- gsub(pattern = "_NT$", 
      #                              replacement = "_CD4", 
      #                              x =  names(results_names), 
      #                              ignore.case = T)
      # 
      # names(results_names) <- gsub(pattern = "_T$", 
      #                              replacement = "_Treg", 
      #                              x =  names(results_names), 
      #                              ignore.case = T)
      
      names(results_names) <- gsub(pattern = "lymphnode", 
                                   replacement = "LN", 
                                   x =  names(results_names), 
                                   ignore.case = T)
      
      names(results_names) <- gsub(pattern = "intraepitheliallayer", 
                                   replacement = "IEL", 
                                   x =  names(results_names), 
                                   ignore.case = T)
      
      names(results_names) <- gsub(pattern = "laminaproprialymphocyte", 
                                   replacement = "LPL", 
                                   x =  names(results_names), 
                                   ignore.case = T)
      
      names(results_names) <- gsub(pattern = "peyerspatches", 
                                   replacement = "PP", 
                                   x =  names(results_names), 
                                   ignore.case = T)
      
      names(results_names) <- gsub(pattern = "_", 
                                   replacement = " ", 
                                   x =  names(results_names), 
                                   ignore.case = T)
    }
    
    ## flow repository IDs ----
    flow_data_ids <- tibble::tibble(
      Dataset = c("Ageing", "Microbiome"), 
      `FlowRepository ID` = c("FR-FCM-Z6L9", "FR-FCM-Z6LT")
    )
    
    ## tSNE info text ----
    tSNSE_info_text <- paste0(
      "The panel included the markers CD103, CD4, CD45, CD62L, CD8a, CD152 (CTLA-4), CD25, CD44, ICOS, CD3, PD-1, CD19, KLRG1, TCR-beta, CD304 (Neuropilin), T-bet, Helios, CD69, NK1.1, ST2, Foxp3, Ki67 and viability. Tregs were gated as viable CD45", tags$sup("+"), 
      " CD3", tags$sup("+"),
      " CD4", tags$sup("+"), 
      " TCR", tags$sup("-"),
      " beta", tags$sup("+"),
      " Foxp3", tags$sup("+"), 
      " CD8", tags$sup("-"),
      " CD19", tags$sup("-"),
      " lymphocytes, and the tSNE analysis was performed on CD103, CTLA-4, CD25, CD44, ICOS, PD-1, KLRG1, Neuropilin, T-bet, Helios, CD69, ST2 and Ki67 using the Cross-Entropy test script in R. The script is ",
      a("available from GitHub.", href = 'https://github.com/AdrianListon/Cross-Entropy-test/tree/main', target = "_blank", .noWS = "outside")
    )
    
    #meta data
    MGI_info <- read.delim(file = "data/combined_MGI_info.txt", 
                           sep = "\t", 
                           stringsAsFactors = F)
    
    #counts
    mean_norm_counts <- readRDS(file = "data/Mean_Normalised_Counts.rds")
    mean_norm_counts$Gene.stable.ID <- gsub(pattern = "\\.[0-9]+$", 
                                            replacement = "", 
                                            x = rownames(mean_norm_counts), 
                                            ignore.case = T)
    #mean_norm_counts2 <- readRDS(file = "Mean_Normalised_Counts.rds")
    #aesthetic formatting of mean_norm_counts
    {
      colnames(mean_norm_counts) <- gsub(pattern = "_NT$", 
                                         replacement = "_Tconv", 
                                         x =  colnames(mean_norm_counts), 
                                         ignore.case = T)
      
      colnames(mean_norm_counts) <- gsub(pattern = "_T$", 
                                         replacement = "_Treg", 
                                         x =  colnames(mean_norm_counts), 
                                         ignore.case = T)
      
      colnames(mean_norm_counts) <- gsub(pattern = "lymphnode", 
                                         replacement = "LN", 
                                         x =  colnames(mean_norm_counts), 
                                         ignore.case = T)
      
      colnames(mean_norm_counts) <- gsub(pattern = "intraepitheliallayer", 
                                         replacement = "IEL", 
                                         x =  colnames(mean_norm_counts), 
                                         ignore.case = T)
      
      colnames(mean_norm_counts) <- gsub(pattern = "laminaproprialymphocyte", 
                                         replacement = "LPL", 
                                         x =  colnames(mean_norm_counts), 
                                         ignore.case = T)
      
      colnames(mean_norm_counts) <- gsub(pattern = "peyerspatches", 
                                         replacement = "PP", 
                                         x =  colnames(mean_norm_counts), 
                                         ignore.case = T)
    }
    
    norm_counts <- DESeq2::counts(object = DESeq_object, normalized = T)
    #aesthetic formatting of norm_counts
    {
      colnames(norm_counts) <- gsub(pattern = "_NT", 
                                    replacement = "_Tconv", 
                                    x =  colnames(norm_counts), 
                                    ignore.case = T)
      
      colnames(norm_counts) <- gsub(pattern = "_T", 
                                    replacement = "_Treg", 
                                    x =  colnames(norm_counts), 
                                    ignore.case = T)
      
      colnames(norm_counts) <- gsub(pattern = "lymphnode", 
                                    replacement = "LN", 
                                    x =  colnames(norm_counts), 
                                    ignore.case = T)
      
      colnames(norm_counts) <- gsub(pattern = "intraepitheliallayer", 
                                    replacement = "IEL", 
                                    x =  colnames(norm_counts), 
                                    ignore.case = T)
      
      colnames(norm_counts) <- gsub(pattern = "laminaproprialymphocyte", 
                                    replacement = "LPL", 
                                    x =  colnames(norm_counts), 
                                    ignore.case = T)
      
      colnames(norm_counts) <- gsub(pattern = "peyerspatches", 
                                    replacement = "PP", 
                                    x =  colnames(norm_counts), 
                                    ignore.case = T)
    }
    
    #PCA, tSNE and UMAP plots
    PCA_means_plots <- readRDS(file = "data/PCA_plots_means.rds")
    PCA_plots <- readRDS(file = "data/PCA_plots.rds")
    tSNE_plots <- readRDS(file = "data/tSNE_plots.rds")
    UMAP_plots <- readRDS(file = "data/UMAP_plots.rds")
    
    #importing min. transcript per cell gene list object
    filtered_genes_lists <- readRDS(file = "data/Filtered_Gene_Lists.rds")
    
    #importing TCR plots
    TCR_plots_list <- list.files(path = "www/TCR_plots/", pattern = "\\.png$")
    #print(TCR_plots_list)
    
    ## aging_tsne_plots_list ----
    aging_tsne_plots_list <- list.files(path = "www/Aging_tsne_plots/", pattern = ".png")
    aging_tsne_tissues <- gsub(pattern = " .*", replacement = "", x = aging_tsne_plots_list, ignore.case = T)
    aging_tsne_tissues <- unique(gsub(pattern = "tsne_plot_", replacement = "", x = aging_tsne_tissues, ignore.case = T))
    
    ## microbiome_tsne_plots_list ----
    microbiome_tsne_plots_list <- list.files(path = "www/Microbiome_tsne_plots/", pattern = ".png")
    microbiome_tsne_tissues <- gsub(pattern = " .*", replacement = "", x = microbiome_tsne_plots_list, ignore.case = T)
    microbiome_tsne_tissues <- unique(gsub(pattern = "tsne_plot_", replacement = "", x = microbiome_tsne_tissues, ignore.case = T))
    
    #importing Pathview plots
    {
      Pathview_plots_list <- list.files(path = "www/Pathview_plots/", 
                                        pattern = "\\.multi\\.png$")
      #print(Pathview_plots_list)
      
      #mapping Pathview plot files to pathway names
      Pathview_mapping <- read.delim(file = "data/output_gage.tab", 
                                     sep = "\t", stringsAsFactors = F)[-1,1:2]
      Pathview_mapping <- 
        Pathview_mapping[Pathview_mapping[,1] %in% 
                           gsub(pattern = "\\..*", 
                                replacement = "", 
                                x = Pathview_plots_list, 
                                ignore.case = T),]
      Pathview_mapping <- Pathview_mapping[order(Pathview_mapping[,1], 
                                                 decreasing = F),]
      
    }
    
  }
  
  #importing functions
  source("R/customFunctions.R")
  
}


#shiny app ui 
{
  # Define UI ----
  ui <- fluidPage(
    
    # using shinyjs to enable or disable certain sections/elements of the app
    shinyjs::useShinyjs(),
    # using shinyalert to enable popup messages
    # shinyalert::useShinyalert(),
    
    # Step-by-step tutorial
    introjsUI(),
    
    # App title ----
    titlePanel("The Tissue Treg Project"),
    
    tabsetPanel(
      
      #Landing Page panel
      {
        tabPanel(
          title = "The Study", 
          value = "TheStudy", 
          h3("High dimensional analysis of tissue-resident regulatory T cells"),
          p("A transcriptional and flow cytometric analysis of tissue-resident regulatory T cells (Tregs) as part of the Tissue Treg Project from the ", 
            a("Liston-Dooley lab.", href='https://www.listonlab.uk/', target = "_blank", .noWS = "outside"),
            " The project is based on Tregs and conventional T cells (Tconv) isolated from blood, spleen, lymph nodes, kidney, liver, lung, pancreas, and gut (Peyer’s patches, LPLs, IELs). Additional tissues covered only in the flow cytometry data sets include adrenals, bone marrow, white adipose tissue (WAT), mesenteric lymph nodes (MLN), brain, skin and thymus. The expression viewer, allowing interactive custom analysis and data download, was developed by Samar Tareen, is ",
            a("open source,", href = 'https://github.com/AdrianListon/ExpressionViewer', 
              target = "_blank", .noWS = "outside"), 
            " and is hosted by the Babraham Institute ",
            a("Bioinformatics Group.", href = 'https://www.bioinformatics.babraham.ac.uk/', 
              target = "_blank", .noWS = "outside"),
            .noWS = c("after-begin", "before-end")
          ),
          hr(),
          column(
            width = 6,
            h4("Flow cytometry data are available for download on FlowRepository using the details below."),
            DT::datatable(
              flow_data_ids,
              rownames = FALSE, 
              options = list(dom = "t")
            )
          ), 
          column(
            width = 4, offset = 2,
            wellPanel(
              h4("Supplementary material for download"),
              
              h5("Bulk RNA-seq analysis"),
              
               withTags(
                ul(
                  li(downloadLink(
                    outputId = "MultiQCBulk", 
                    label = "QC and pre-processing report"
                  )),
                  li(downloadLink(
                    outputId = "GlobalBulk", 
                    label = "PCA, t-SNE and UMAP plots"
                  )),
                  li(downloadLink(
                    outputId = "DiffExprsBulk", 
                    label = paste0("Counts and differential ", 
                                   "expression results")
                  )),
                  li(downloadLink(
                    outputId = "GSEABulk", 
                    label = "Gene set enrichment results"
                  ))
                )
              )
            )
          ),
          
          br(),
          br(),
          br(),
          br(),
          br(),
          
          ## footer text ----
          fluidRow(
            h6(
              style="padding:10px;", 
              HTML(
                "Funded by the European Union (ERC, TissueTreg, 681373)
                 and Wellcome Trust (Brain CD4 T cells and their influence 
                 over microglial homeostasis, 222442/A/21/Z). <br>
                 Copyright 2023, Liston-Dooley Lab, 
                 Babraham Institute and University of Cambridge."
              )
            )
          )
        )
      },
      
      #Global Analysis panel
      {
        tabPanel(title = "Global Analysis", value = "GlobalAnalysis",
                 
                 h3("Global analysis of tissue Treg transcriptome across different tissue sources"),
                 
                 p(paste0("For the bulk RNA-seq analysis, Tregs and Tconv were ", 
                          "isolated from blood, and perfused spleen, lymph ", 
                          "nodes, kidney, liver, lung, pancreas, and gut ", 
                          "(Peyer’s patches, LPLs, IELs) samples. Samples were ", 
                          "prepared from Foxp3Thy1.1 reporter mice, allowing ", 
                          "purification based on Foxp3 expression. mRNA was ", 
                          "sequencing using the QuantSeq 3’mRNA-Seq Library ", 
                          "Prep Kit for Illumina and the QuantSeg data analysis ", 
                          "workflow. The relationship between samples at the ", 
                          "global level is visualised through PCA (using all ", 
                          "samples or tissue means), tSNE and UMAP projections.")),
                 
                 # p(paste0("For the single cell RNA-seq analysis, Tregs were ", 
                 #          "flow sorted from Foxp3Thy1.1 reporter mice, ", 
                 #          "pre-injected with anti-CD45 via intravenous ", 
                 #          "delivery. Cells were purified from blood, kidney, ", 
                 #          "liver, pancreas and LPLs from the gut. Purified ", 
                 #          "cells were CD4+Foxp3Thy1.1+ as well as negative for ", 
                 #          "intravenous CD45 labelling and for the exclusion ", 
                 #          "markers CD19, CD11b, CD8 and F4/80. Cells were ", 
                 #          "labelled with Hastag TotalSeq and loaded onto the ", 
                 #          "10x Chromium Controller with sequencing performed ", 
                 #          "on Illumina HiSeq. Single cell data is visually ", 
                 #          "presented as tSNE and UMAP projections, with each ", 
                 #          "tissue origin labelled a different colour.")),
                 
                 p(paste0("Plots may be downloaded by right clicking the ",
                          "plot and selecting 'Save image as...'")),
                 
                 #br(),
                 hr(),
                 
                 h4("Bulk RNA-seq data"),
                 
                 fluidRow(
                   
                   column(width = 12,
                          
                          selectizeInput(inputId = "TriplePlot", 
                                         label = "Select Plot Type", 
                                         # choices = as.list(results_names),
                                         #choices = c("", "PCA", "t-SNE", "UMAP"),
                                         choices = c("PCA - tissue means", 
                                                     "PCA", "t-SNE", "UMAP"),
                                         selected = "PCA - tissue means", 
                                         #selectize = T, 
                                         multiple = F, 
                                         options = list(
                                           placeholder = "Click to select")
                          ),
                   )
                 ),
                 
                 fluidRow(
                   
                   column(width = 4,
                          plotOutput(outputId = "PCAplot1") %>% 
                            withSpinner(type = 3, size = 0.5, 
                                        color.background = "#FFFFFF")
                   ),
                   column(width = 4,
                          plotOutput(outputId = "PCAplot2") %>% 
                            withSpinner(type = 3, size = 0.5, 
                                        color.background = "#FFFFFF")
                   ),
                   column(width = 4,
                          plotOutput(outputId = "PCAplot3") %>% 
                            withSpinner(type = 3, size = 0.5, 
                                        color.background = "#FFFFFF")
                   )
                 ),
                 
                 hr(),
                 
                 # h4("Single cell RNA-seq data"),
                 # 
                 # p("Will be available soon")
                 
        )
      },
      
      # Bulk RNA-seq Expression panel ----
      {
        tabPanel(title = "Bulk RNA-seq Expression", 
                 value = "BulkDiffExprs",
                 
                 br(),              ### browser button ----
                 actionButton("browser", "browser"),
                 
                 fluidRow(
                   
                   column(width = 7, 
                          h3("Generate contrasts"),
                          
                          "An interactive differential analysis tool based ", 
                          "on the bulk RNAseq dataset. Pick any population ", 
                          "or group of populations to run a differential ", 
                          "expression analysis by comparison to any other ", 
                          "population or group of populations. Both Tregs ", 
                          "and Tconv are available, as purified from the ", 
                          "blood, spleen, lymph nodes, kidney, liver, lung, ", 
                          "pancreas, Peyer’s Patches, and intraepithelial ", 
                          "cells (IEL) and lamina propria lymphocytes (LPL) ", 
                          "from the intestines. The differential expression ", 
                          "analysis is performed using DESeq2, using a custom ", 
                          "script available on ", 
                          a("GitHub.", href= 'https://github.com/AdrianListon/TissueTregs', 
                            target = "_blank"),
                          p(),
                          "For a quick walkthrough on doing the analysis, ", 
                          "please click the yellow button with the question ", 
                          "mark. Download options for visualization and ", 
                          "data tables become available after running a contrast.",
                          
                   ),
                   
                   column(width = 5, 
                          br(),
                          br(),
                          # actionButton(inputId = "TutorialButton", 
                          #              label = "Start tutorial"),
                          actionBttn(inputId = "TutorialButton", 
                                     label = "Start tutorial", 
                                     style = "material-circle", 
                                     color = "warning", 
                                     icon = icon("question")),
                   )
                 ),
                 
                 hr(),
                 
                 fluidRow(
                   
                   column(width = 2, 
                          
                          wellPanel(
                            
                            fluidRow(
                              
                              column(width = 9, 
                                     h4("Available Contrasts"),
                              ),
                              
                              column(width = 3,
                                     
                                     dropdownButton(inputId = "mydropdown",
                                                    circle = TRUE,
                                                    #status = "info",
                                                    status = "secondary",
                                                    #icon = icon("question"),
                                                    icon = icon("info"),
                                                    right = F, up = F,
                                                    width = "500px",
                                                    tooltip = tooltipOptions(
                                                      title = "Instructions"),
                                                    
                                                    h4("To Begin:"),
                                                    
                                                    h5(paste0("Start by ", 
                                                              "selecting tissues",
                                                              " for the contrast", 
                                                              " groups in this ", 
                                                              "panel. This will", 
                                                              " enable the ", 
                                                              "'Generate ", 
                                                              "Contrasts' button", 
                                                              " at the bottom. ", 
                                                              "A tissue cannot ", 
                                                              "be in both ", 
                                                              "groups ", 
                                                              "simultaneously.")),
                                                    
                                                    h5(paste0("The drop down ",
                                                              "list and 'Add' ",
                                                              "button below ",
                                                              "each of the two ", 
                                                              "selections can ", 
                                                              "be used to add ", 
                                                              "relevant tissues", 
                                                              ", or to ", 
                                                              "clear the current", 
                                                              " selection. The ", 
                                                              "buttons will not ", 
                                                              "duplicate ", 
                                                              "tissues that are ", 
                                                              "already selected."
                                                    )),
                                                    
                                                    h5(paste0("Clicking the ", 
                                                              "'Generate ", 
                                                              "Contrasts' button", 
                                                              " will perform ", 
                                                              "the differential ", 
                                                              "expression ", 
                                                              "analysis between ", 
                                                              "the two contrast ", 
                                                              "groups, creating ",
                                                              "a table of ", 
                                                              "results and a ", 
                                                              "volcano plot on ", 
                                                              "right."
                                                    )),
                                                    
                                                    h5(paste0("Clicking the ", 
                                                              "'Show Advanced Options' ", 
                                                              "checkbox reveals three ", 
                                                              "additional options to ", 
                                                              "adjust the analysis and ", 
                                                              "visualisation. The ", 
                                                              "'Minimum Transcripts per ", 
                                                              "Cell' drop down list ", 
                                                              "gives options to filter ", 
                                                              "results to remove low ", 
                                                              "expression genes if ", 
                                                              "needed. The 'log2 Fold ", 
                                                              "Change' and 'Adjusted ", 
                                                              "p.value' fields are the ", 
                                                              "significance criteria ", 
                                                              "used in the volcano plot ", 
                                                              "and can be adjusted as ", 
                                                              "needed."
                                                    )),
                                                    
                                                    hr(),
                                                    h4("Interactivity:"),
                                                    
                                                    h5(paste0("The table and ", 
                                                              "volcano plot are ", 
                                                              "interactive. Up ", 
                                                              "to 50 points/rows", 
                                                              " can be selected ", 
                                                              "to mark them in ", 
                                                              "the volcano plot ", 
                                                              "and to generate ", 
                                                              "a plot of ", 
                                                              "normalised counts", 
                                                              " of the selected ", 
                                                              "genes across all ", 
                                                              "samples.")),
                                                    
                                                    h5(paste0("Alternatively, ",
                                                              "points can also", 
                                                              " be selected by ", 
                                                              "clicking on them", 
                                                              " in the volcano", 
                                                              " plot to mark ", 
                                                              "them in the plot", 
                                                              " and also select", 
                                                              " the respective", 
                                                              " rows in the ", 
                                                              "table. Clicking", 
                                                              " a selected ", 
                                                              "point again will", 
                                                              " deselect it in ", 
                                                              "the plot and ", 
                                                              "table.")),
                                                    
                                                    h5(paste0("The selection of ", 
                                                              "points/rows can ", 
                                                              "be reset at any ", 
                                                              "time using the ", 
                                                              "'Reset Row ", 
                                                              "Selection' button", 
                                                              " above the table." 
                                                    )),
                                                    
                                                    h5(paste0("The plots can be ", 
                                                              "zoomed in by ", 
                                                              "dragging a ", 
                                                              "selection box ", 
                                                              "over the relevant", 
                                                              " area. Reset the ", 
                                                              "zoom by either ", 
                                                              "double clicking ", 
                                                              "anywhere in the ", 
                                                              "plot or by ", 
                                                              "clicking the home", 
                                                              " icon on the top ", 
                                                              "right of the ", 
                                                              "plot.")),
                                                    
                                                    h5(paste0("The table can be ", 
                                                              "filtered through ", 
                                                              "the fields under ", 
                                                              "the column ", 
                                                              "headers. These ", 
                                                              "can also be used ", 
                                                              "to search the ", 
                                                              "column of ", 
                                                              "interest.")),
                                                    
                                                    h5(paste0("Additionaly, the ", 
                                                              "table can be ", 
                                                              "sorted by ", 
                                                              "multiple columns ", 
                                                              "by shift clicking", 
                                                              " the column ", 
                                                              "headers.")),
                                                    
                                                    hr(),
                                                    h4("Download/Export:"),
                                                    
                                                    h5(paste0("The results table", 
                                                              " can be ", 
                                                              "downloaded as ", 
                                                              "tab delimited ", 
                                                              "files from the ", 
                                                              "three download ", 
                                                              "buttons above the", 
                                                              " table.")),
                                                    
                                                    h5(paste0("The plots can be ", 
                                                              "downloaded via ", 
                                                              "the camera icon ", 
                                                              "on the top right ", 
                                                              "of each plot.", 
                                                              "Downloading a ", 
                                                              "plot while zoomed", 
                                                              " in will download", 
                                                              " the zoomed in ", 
                                                              "view.")),
                                                    
                                                    h5(paste0("The 'Show Plot ", 
                                                              "Export Options' ", 
                                                              "that becomes ", 
                                                              "available above ", 
                                                              "the volcano plots", 
                                                              " reveals options ", 
                                                              "to modify the ", 
                                                              "export type and ", 
                                                              "size of the ", 
                                                              "plot. The ", 
                                                              "format of the ", 
                                                              "image (SVG/PNG/", 
                                                              "JPEG) can be ", 
                                                              "selected from ", 
                                                              "the radio button ", 
                                                              "that becomes ", 
                                                              "available below ",  
                                                              "the 'Generate ", 
                                                              "Conrasts' button.", 
                                                              " The size of the ", 
                                                              "plot can be ", 
                                                              "adjusted by ", 
                                                              "specifying the ", 
                                                              "height and width ", 
                                                              "in pixels."
                                                    ))
                                     )
                              ),
                            ),
                            #h4("Available Contrasts"),
                            
                            selectizeInput(inputId = "Group1", 
                                           label = "Contrast Group 1", 
                                           # choices = as.list(results_names),
                                           choices = list(
                                             Tregs = results_names[
                                               grepl(pattern = "cellTypeT",
                                                     #grepl(pattern = "_Treg",
                                                     x = results_names,
                                                     ignore.case = T)],
                                             Tconv = results_names[
                                               grepl(pattern = "cellTypeNT",
                                                     #grepl(pattern = "_Tconv",
                                                     x = results_names,
                                                     ignore.case = T)]
                                           ),
                                           selected = NULL, 
                                           #selectize = T, 
                                           multiple = T, 
                                           options = list(
                                             placeholder = paste0("Click here ", 
                                                                  "to manually",
                                                                  " select ", 
                                                                  "...")
                                           )
                                           #)
                            ),
                            
                            #fluidRow(
                            #h5("Select available:"),
                            #h5("Quick add available:"),
                            # actionButton(inputId = "Group1SelectTconv", 
                            #              label = "Tconv"),
                            # actionButton(inputId = "Group1SelectTregs", 
                            #              label = "Tregs"),
                            selectizeInput(inputId = "Group1Selector",
                                           label = "Add available:",
                                           choice = list(
                                             Treg = c("All Tregs", 
                                                      "Blood Tregs", 
                                                      "Lymphoid Tregs", 
                                                      "Non-lymphoid Tregs", 
                                                      "Gut-associated Tregs"),
                                             Tconv = c("All Tconv", 
                                                       "Blood Tconv", 
                                                       "Lymphoid Tconv", 
                                                       "Non-lymphoid Tconv", 
                                                       "Gut-associated Tconv")
                                           ), 
                                           selected = "All Tregs", 
                                           multiple = F, 
                                           options = list(
                                             placeholder = "Click to select")
                            ),
                            actionButton(inputId = "Group1Add", 
                                         label = "Add"),
                            actionButton(inputId = "Group1Clear", 
                                         label = "Clear Selection"),
                            #),
                            
                            #br(),
                            #br(),
                            hr(),
                            #br(),
                            
                            #output of group 1 selection for debugging -- remove/
                            # hide when deploying
                            #verbatimTextOutput(outputId = "value_Group1"),
                            
                            #p("compared against"),
                            
                            #selection for group 2 (denominator) of the contrasts
                            #initially disabled to allow 
                            #disabled(
                            selectizeInput(inputId = "Group2", 
                                           label = "Contrast Group 2", 
                                           choices = list(
                                             Tregs = results_names[
                                               grepl(pattern = "cellTypeT",
                                                     x = results_names,
                                                     ignore.case = T)],
                                             Tconv = results_names[
                                               grepl(pattern = "cellTypeNT",
                                                     x = results_names,
                                                     ignore.case = T)]
                                           ),
                                           selected = NULL, 
                                           #selectize = T, 
                                           multiple = T, 
                                           options = list(
                                             placeholder = paste0("Click here ", 
                                                                  "to manually",
                                                                  " select ", 
                                                                  "...")
                                           )
                            ),
                            
                            #fluidRow(
                            #h5("Select available:"),
                            #h5("Quick add available:"),
                            # actionButton(inputId = "Group2SelectTconv", 
                            #              label = "Tconv"),
                            # actionButton(inputId = "Group2SelectTregs", 
                            #              label = "Tregs"),
                            selectizeInput(inputId = "Group2Selector",
                                           label = "Add available:",
                                           choice = list(
                                             Treg = c("All Tregs", 
                                                      "Blood Tregs", 
                                                      "Lymphoid Tregs", 
                                                      "Non-lymphoid Tregs", 
                                                      "Gut-associated Tregs"),
                                             Tconv = c("All Tconv", 
                                                       "Blood Tconv", 
                                                       "Lymphoid Tconv", 
                                                       "Non-lymphoid Tconv", 
                                                       "Gut-associated Tconv")
                                           ), 
                                           selected = "All Tconv", 
                                           multiple = F, 
                                           options = list(
                                             placeholder = "Click to select")
                            ),
                            actionButton(inputId = "Group2Add", 
                                         label = "Add"),
                            actionButton(inputId = "Group2Clear", 
                                         label = "Clear Selection"),
                            #),
                            
                            hr(),
                            
                            checkboxInput(inputId = "AddOps", 
                                          label = "Show Additional Options"),
                            includeScript(path = "js/validate-number.js"),
                            conditionalPanel(
                              condition = "input.AddOps",
                              selectizeInput(inputId = "TPCFilt", 
                                             label = "Minimum Transcripts per Cell",
                                             choice = c("No filter", "0.01", 
                                                        "0.1", "1"), 
                                             selected = "0.01", 
                                             multiple = F),
                              h5("Volcano Plot Marking Criteria"),
                              numericInput(inputId = "log2FCFilt",
                                           label = HTML(paste0("log2 Fold ", 
                                                               "Change (+/-)", 
                                                               "<br/>(Range: ", 
                                                               "0 - 20)")),
                                           value = 2,
                                           min = 0,
                                           max = 20,
                                           step = 0.1),
                              numericInput(inputId = "AdjpFilt",
                                           label = "Adjusted p.value", 
                                           value = 0.01, 
                                           min = 0, 
                                           max = 1, 
                                           step = 0.01),
                            ),
                            hr(),
                            #output of group 2 selection for debugging -- remove/
                            # hide when deploying
                            #verbatimTextOutput(outputId = "value_Group2"),
                            
                            #action button for starting the DESeq2 results() 
                            # function for generating contrasts
                            actionButton(inputId = "GenContrasts",
                                         label = "Generate Contrasts"
                            ),
                            
                          ),
                          
                          # verbatimTextOutput(outputId = "res_error"),
                          # 
                          # h6("Debug console-like thingy"),
                          # verbatimTextOutput(outputId = "debug_console")
                          
                   ),
                   
                   column(width = 4, 
                          
                          conditionalPanel(
                            condition = "output.show_panels",
                            
                            #h5(HTML("<br/>")),
                            br(),
                            fluidRow(
                              column(width = 10,
                                     p(paste0("For an overview of the results, ", 
                                              "features and interactability, please ", 
                                              "click the yellow button to the right."))
                              ),
                              column(width = 2,
                                     # actionButton(inputId = "TutorialButton2", 
                                     #              label = "Start tutorial"),
                                     actionBttn(inputId = "TutorialButton2", 
                                                label = "Start tutorial", 
                                                style = "material-circle", 
                                                color = "warning", 
                                                icon = icon("question")),
                              )
                            ),
                            
                            hr(),
                            
                            checkboxInput(inputId = "PlotSaveOps", 
                                          label = "Show Plot Export Options"), 
                            
                            conditionalPanel(
                              condition = "input.PlotSaveOps",
                              
                              fluidRow(
                                column(width = 3,#)
                                       radioGroupButtons(inputId = "SavePlotFormat", 
                                                         label = HTML(paste0("Plot", 
                                                                             " Export ", 
                                                                             "Format")), 
                                                         choices = c("SVG", "PNG", "JPEG"), 
                                                         #SVG, PNG, JPEG and WEBP 
                                                         #np native PDF support
                                                         selected = "SVG", 
                                                         direction = "vertical", #"horizontal", 
                                                         checkIcon = list(
                                                           # yes = icon("ok", 
                                                           #            lib = "glyphicon")
                                                           yes = icon(name = "file-image", 
                                                                      lib = "font-awesome")
                                                         ),
                                                         size = "normal" # xs/sm/normal/lg
                                       ),
                                ),
                                column(width = 7,#)
                                       includeScript(path = "js/validate-number.js"),
                                       
                                       numericInputIcon(inputId = "WidthPixels", 
                                                        label = HTML(paste0("Plot ", 
                                                                            "Size in ",  
                                                                            " Pixels ", #"<br/>", 
                                                                            "(Range: 100 ", 
                                                                            "- 10,000) ", 
                                                                            #"<br/><br/>", 
                                                                            "<br/>", 
                                                                            "Width")), #in Pixels",
                                                        value = 880,
                                                        min = 100,
                                                        max = 10000,
                                                        step = 10,
                                                        width = "85%",
                                                        icon = list(
                                                          icon(name = "ruler-horizontal",
                                                               lib = "font-awesome"),
                                                          "px")
                                       ),
                                       
                                       numericInputIcon(inputId = "HeightPixels", 
                                                        label = "Height", #in Pixels",
                                                        value = 562,
                                                        min = 100,
                                                        max = 10000,
                                                        step = 10,
                                                        width = "85%",
                                                        icon = list(
                                                          icon(name = "ruler-vertical",
                                                               lib = "font-awesome"),
                                                          "px")
                                       )
                                )
                              )
                            ),
                            
                            hr(),
                            
                            h4("Differential Expression Volcano Plot"),
                            # h5(paste0("Click on data points to select/deselect ",
                            #           "them in the table and show/hide their ",
                            #           "normalised counts.")),
                            
                            plotlyOutput(outputId = "res_volplot") %>% 
                              withSpinner(type = 3, size = 0.5, 
                                          color.background = "#FFFFFF"),
                          ),
                          
                          conditionalPanel(
                            condition = paste0(
                              "typeof input.res_table_rows_selected  !== ", 
                              "'undefined' && ", 
                              "input.res_table_rows_selected.length > 0"
                            ),
                            hr(),
                            h4("Normalised Expression Plot"),
                            
                            plotlyOutput(outputId = "res_colplot") %>% 
                              withSpinner(type = 3, size = 0.5, 
                                          color.background = "#FFFFFF")
                          )
                          
                   ),
                   
                   column(width = 5,
                          
                          fluidRow(
                            conditionalPanel(
                              condition = "output.show_panels",
                              
                              #fluidRow(
                              #column(width = 11,
                              h4(HTML("<br/>Results Table")),
                              h5("Download tab delimited file:"),
                              downloadButton(
                                outputId = "DownloadTableAll", 
                                label = "Complete Results"
                              ),
                              downloadButton(
                                outputId = "DownloadTableFilt", 
                                label = "Filtered Table"
                              ),
                              downloadButton(
                                outputId = "DownloadTableSelected", 
                                label = "Selected Rows Only"
                              )
                              #)
                              #),
                            ),
                          ),
                          #br(),
                          conditionalPanel(
                            condition = "output.show_panels",
                            fluidRow(
                              column(width = 12,
                                     hr(),
                                     #h5("Row Selection Options:"),
                                     actionButton(inputId = "ResetTable", 
                                                  label = "Reset Row Selection"
                                     )
                              )
                            ),
                            br(),
                          ),
                          
                          # conditionalPanel(
                          #   condition = paste0(
                          #     "typeof output.res_table  !== 'undefined'"
                          #     ),
                          #shinyjs::hidden(
                          fluidRow(
                            #shinyjs::hidden(
                            DT::dataTableOutput(outputId = "res_table") #%>% 
                            # withSpinner(type = 3, size = 0.5, 
                            #             color.background = "#FFFFFF")
                            #)
                          )
                   )
                 )
        )
      },
      
      # Pathway analysis panel ----
      {
        tabPanel(title = "Pathway Analysis", value = "PathwayAnalysis", 
                 
                 h3("Visualise enriched KEGG pathways"),
                 
                 p(paste0("Tregs and Tconv were isolated from blood, and ", 
                          "perfused spleen, lymph nodes, kidney, liver, lung, ", 
                          "pancreas, and gut (Peyer’s patches, LPLs, IELs) ", 
                          "samples. Samples were prepared from Foxp3Thy1.1 ", 
                          "reporter mice, allowing purification based on Foxp3 ", 
                          "expression. Bulk mRNA was sequencing using the ", 
                          "QuantSeq 3’mRNA-Seq Library Prep Kit for Illumina ", 
                          "and the QuantSeg data analysis workflow. Differential ", 
                          "analysis was performed between Tregs in the blood, ", 
                          "lymphoid tissues (spleen, lymph nodes), non-lymphoid ", 
                          "organs (kidney, liver, lung, pancreas) and ", 
                          "gut-associated tissues (Peyer’s patches, LPLs, IELs). ", 
                          "KEGG pathways enriched for differential expression ", 
                          "are visualised using GAGE and Pathview, with each ", 
                          "pairwise differential expression set visualised as ", 
                          "log2 fold-change for any gene within the pathway with ", 
                          "differential expression.")),
                 
                 hr(),
                 
                 #vr(),
                 
                 column(width = 9, 
                        
                        fluidRow(
                          
                          #h4("Visualise Enriched KEGG Pathways"), 
                          
                          column(width = 6, 
                                 #br(),
                                 br(),
                                 
                                 selectizeInput(inputId = "PathviewPlotSelect", 
                                                label = "Select KEGG Pathway", 
                                                # choices = as.list(results_names),
                                                #choices = c("", "PCA", "t-SNE", "UMAP"),
                                                choices = paste(Pathview_mapping[,1], 
                                                                Pathview_mapping[,2], 
                                                                sep = " - "),
                                                selected = paste(Pathview_mapping[,1], 
                                                                 Pathview_mapping[,2], 
                                                                 sep = " - ")[1], 
                                                #selectize = T, 
                                                multiple = F, 
                                                options = list(
                                                  placeholder = "Click to select"), 
                                                width = "600px"
                                 ),
                                 
                                 downloadButton(
                                   outputId = "DownloadPathviewSelected", 
                                   label = "Selected Pathway"
                                 ),
                                 
                                 downloadButton(
                                   outputId = "DownloadPathviewAll", 
                                   label = "All Pathways"
                                 )
                                 
                          ),
                          
                          column(width = 3, 
                                 
                                 imageOutput(outputId = "PathviewPlotLegend", 
                                             width = "100%", 
                                             height = "100%") %>% 
                                   withSpinner(type = 3, size = 0.5, 
                                               color.background = "#FFFFFF"),
                                 br()
                          ),
                          
                          column(width = 3, 
                                 
                                 br(),
                                 imageOutput(outputId = "PathviewPlotLegend2", 
                                             width = "100%", 
                                             height = "100%") %>% 
                                   withSpinner(type = 3, size = 0.5, 
                                               color.background = "#FFFFFF")
                          )#,
                        ),
                        
                        imageOutput(outputId = "PathviewPlotImage", 
                                    width = 1342*0.35, 
                                    height = 1451*0.35) %>% 
                          withSpinner(type = 3, size = 0.5, 
                                      color.background = "#FFFFFF")
                        
                 )
        )
        
      },
      
      #Single Cell Expression panel
      # {
      #   tabPanel(title = "Single Cell Expression", value = "SingleCellExprs", 
      #            
      #            )
      # },
      
      # TCR panel ----
      {
        tabPanel(title = "TCR Repertoire", value = "TCR", 
                 
                 h3("TCR repertoire analysis"),
                 
                 p(paste0("Tregs were flow sorted from Foxp3Thy1.1 reporter ", 
                          "mice, pre-injected with anti-CD45 via intravenous ", 
                          "delivery. Cells were purified from blood, kidney, ", 
                          "liver, pancreas and LPLs from the gut. Purified ", 
                          "cells were CD4+Foxp3Thy1.1+ as well as negative for ", 
                          "intravenous CD45 labelling and for the exclusion ", 
                          "markers CD19, CD11b, CD8 and F4/80. Cells were ", 
                          "labelled with Hastag TotalSeq and loaded onto the ", 
                          "10x Chromium Controller with sequencing performed on ", 
                          "Illumina HiSeq. Data is visualised for each cell ", 
                          "where TCRbeta CDR3 expression was observed using ", 
                          "Chord diagrams. Chord diagrams represent amino ", 
                          "acid-level unique clonotypes on the circular axis, ", 
                          "with shared clonotypes represented on the radial ", 
                          "axis as shared across different tissues in the same ", 
                          "mouse (green), shared across different tissues in ", 
                          "different mice (black) and shared between the same ", 
                          "tissue in different mice (pink).")),
                 
                 #br(),
                 hr(),
                 
                 fluidRow(
                   
                   column(width = 4,
                          
                          selectizeInput(inputId = "TCRChordSelect", 
                                         label = "Select Tissue", 
                                         # choices = as.list(results_names),
                                         #choices = c("", "PCA", "t-SNE", "UMAP"),
                                         choices = gsub(pattern = "_.*", 
                                                        replacement = "", 
                                                        x = TCR_plots_list, 
                                                        ignore.case = T),
                                         selected = "Blood", 
                                         #selectize = T, 
                                         multiple = F, 
                                         options = list(
                                           placeholder = "Click to select")
                          ),
                          
                          downloadButton(
                            outputId = "DownloadTCRSelected", 
                            label = "Selected Tissue"
                          ),
                          
                          downloadButton(
                            outputId = "DownloadTCRAll", 
                            label = "All Tissues"
                          )
                          
                   ),
                   
                   column(width = 8,
                          
                          imageOutput(outputId = "TCRChordImage", 
                                      width = 4344*0.15, 
                                      height = 4601*0.15) %>% 
                            withSpinner(type = 3, size = 0.5, 
                                        color.background = "#FFFFFF")
                          
                   )
                   
                 )
        )
      },
      
      # Aging t-SNE panel ----
      {
        tabPanel(
          title = "Aging t-SNE", 
          value = "AgingTSNE", 
          h3("t-SNE plots across 6 time-points"),
          p(
            HTML(
              paste0("Leukocytes were isolated from the perfused tissues of mice aged 8 to 100 weeks, stained for flow cytometry and acquired on a BD FACSymphony. ", tSNSE_info_text)
              # paste0(
              #   "Leukocytes were isolated from the perfused tissues of mice aged 8 to 100 weeks, stained for flow cytometry and acquired on a BD FACSymphony. The panel included the markers CD103, CD4, CD45, CD62L, CD8a, CD152 (CTLA-4), CD25, CD44, ICOS, CD3, PD-1, CD19, KLRG1, TCR-beta, CD304 (Neuropilin), T-bet, Helios, CD69, NK1.1, ST2, Foxp3, Ki67 and viability. Tregs were gated as viable CD45", 
              # tags$sup("+"), 
              # " CD3",
              # tags$sup("+"),
              # " CD4",
              # tags$sup("+"), 
              # " TCR",
              # tags$sup("-"),
              # " beta",
              # tags$sup("+"),
              # " Foxp3",
              # tags$sup("+"), 
              # " CD8",
              # tags$sup("-"),
              # " CD19",
              # tags$sup("-"),
              # " lymphocytes, and the tSNE analysis was performed on CD103, CTLA-4, CD25, CD44, ICOS, PD-1, KLRG1, Neuropilin, T-bet, Helios, CD69, ST2 and Ki67 using the Cross-Entropy test script in R. The script is ",
              # a("available from GitHub,", href = 'https://github.com/AdrianListon/Cross-Entropy-test/tree/main', target = "_blank", .noWS = "outside")
              #)
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 2,
              selectizeInput(
                inputId = "agingTSNEtissue", 
                label = "Select Tissue", 
                choices = aging_tsne_tissues,
                #selectize = T, 
                multiple = FALSE, 
                options = list(
                  placeholder = "Click to select")
              ),
            ),
            column(
              width = 2, 
              offset = 1,
              br(),
              downloadButton(
                outputId = "DownloadAgingTSNEtissue", 
                label = "Selected Tissue"
              )
            ),
            column(
              width = 1,
              br(),
              downloadButton(
                outputId = "DownloadAgingTSNEAll", 
                label = "All Tissues"
              )
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 4,
              textOutput(outputId = "agingTSNETitle1"),
              imageOutput(outputId = "agingTSNEImage1", width = "100%", height = "100%")
            ),
            column(
              width = 4,
              textOutput(outputId = "agingTSNETitle2"),
              imageOutput(outputId = "agingTSNEImage2", width = "100%", height = "100%")
            ),
            column(
              width = 4,
              textOutput(outputId = "agingTSNETitle3"),
              imageOutput(outputId = "agingTSNEImage3", width = "100%", height = "100%") 
            )
          ),
          br(),
          fluidRow( 
            column(
              width = 4,
              textOutput(outputId = "agingTSNETitle4"),
              imageOutput(outputId = "agingTSNEImage4", width = "100%", height = "100%")
            ),
            column(
              width = 4,
              textOutput(outputId = "agingTSNETitle5"),
              imageOutput(outputId = "agingTSNEImage5", width = "100%", height = "100%")
            ),
            column(
              width = 4,
              textOutput(outputId = "agingTSNETitle6"),
              imageOutput(outputId = "agingTSNEImage6", width = "100%", height = "100%")
            )
          ) 
        )
      },
      
      # Microbiome t-SNE panel ----
      {
        tabPanel(
          title = "Microbiome t-SNE", 
          value = "MicrobiomeTSNE", 
          h3("t-SNE plots across 3 conditions"),
          p(
            HTML(
              paste0(
            "The impact of the microbiome on Treg phenotype was assessed by high parameter flow cytometry. SPF mice were compared to gnotobiotic (germ-free) and wilded (cohoused) microbiome mice. For the microbiome enrichment, pet store female mice were wild-exposed prior to cohousing with SPF C57BL/6J mice. Leukocytes were isolated from the perfused tissues of mice, stained for flow cytometry and acquired on a BD FACSymphony. ",
            tSNSE_info_text
            ))
          ),
          hr(),
          fluidRow(
            column(
              width = 2,
              selectizeInput(
                inputId = "microbiomeTSNEtissue", 
                label = "Select Tissue", 
                choices = microbiome_tsne_tissues,
                #selectize = T, 
                multiple = FALSE, 
                options = list(
                  placeholder = "Click to select")
              ),
            ),
            column(
              width = 2, 
              offset = 1,
              br(),
              downloadButton(
                outputId = "DownloadMicrobiomeTSNEtissue", 
                label = "Selected Tissue"
              )
            ),
            column(
              width = 1,
              br(),
              downloadButton(
                outputId = "DownloadMicrobiomeTSNEAll", 
                label = "All Tissues"
              )#,
              #actionButton("browser", "browser")
            )
          ),
          hr(),
          fluidRow(
            column(
              width = 4,
              textOutput(outputId = "microbiomeTSNETitle1"),
              imageOutput(outputId = "microbiomeTSNEImage1", width = "100%", height = "100%")
            ),
            column(
              width = 4,
              textOutput(outputId = "microbiomeTSNETitle2"),
              imageOutput(outputId = "microbiomeTSNEImage2", width = "100%", height = "100%")
            ),
            column(
              width = 4,
              textOutput(outputId = "microbiomeTSNETitle3"),
              imageOutput(outputId = "microbiomeTSNEImage3", width = "100%", height = "100%") 
            )
          )
        )
      }
    )
  )
  ## ui fluid page ends at bracket immediately above
  
}


#shiny app server
{
  
  # Define server logic ----
  server <- function(input, output, session) {
    
    observeEvent(input$browser, browser())
    
    #reading inputs for debugging
    #output$value_Group1 <- renderPrint({ (input$Group1) })
    #output$value_Group2 <- renderPrint({ (input$Group2) })
    
    #output$value_padj <- renderPrint({ (input$Padj_cutoff) })
    
    #debug console output
    # output$debug_console <- renderPrint(results_names[
    #   grepl(pattern = "_NT",
    #         x = results_names,
    #         ignore.case = T)])
    # output$debug_console <- renderPrint(class(input$Group1))
    
    #external variable handles so that they do not get destroyed once relevant
    #gui updates and rendering are performed
    global_data <- reactiveValues(result = NULL, 
                                  base_traces = 0, 
                                  #plotly traces are 0 indexed so count is ct-1
                                  trace_count = 0, 
                                  #trace_count = base_traces-1, 
                                  dt = NULL,
                                  #volplot_selected = NULL#,
                                  #volplot = NULL,
                                  #colplot = NULL
                                  left_log = 0, 
                                  right_log = 0,
                                  first_run_done = F
    )
    
    #live updating UI elements
    observe({
      
      updateSelectInput(session = session,
                        inputId = "Group1",
                        label = "Contrast Group 1",
                        # choices = as.list(results_names[
                        #   !(results_names %in% input$Group2)
                        #   ]),
                        choices = list(
                          Tregs = results_names[
                            grepl(pattern = "cellTypeT",
                                  x = results_names,
                                  ignore.case = T) & !(results_names %in%
                                                         input$Group2)],
                          Tconv = results_names[
                            grepl(pattern = "cellTypeNT",
                                  x = results_names,
                                  ignore.case = T) & !(results_names %in%
                                                         input$Group2)]
                        ),
                        selected = input$Group1)
      
      updateSelectInput(session = session,
                        inputId = "Group2",
                        label = "Contrast Group 2",
                        # choices = as.list(results_names[
                        #   !(results_names %in% input$Group1)
                        #   ]),
                        choices = list(
                          Tregs = results_names[
                            grepl(pattern = "cellTypeT",
                                  x = results_names,
                                  ignore.case = T) & !(results_names %in%
                                                         input$Group1)],
                          Tconv = results_names[
                            grepl(pattern = "cellTypeNT",
                                  x = results_names,
                                  ignore.case = T) & !(results_names %in%
                                                         input$Group1)]
                        ),
                        selected = input$Group2)
      
      updateNumericInput(session = session, 
                         inputId = "HeightPixels", 
                         # label = tags$div(HTML(paste0('<i class="fa fa-ruler-', 
                         #                              'vertical"></i> Plot ', 
                         #                              'Height (pixels)')
                         #                       )
                         #                  ), 
                         value = ifelse(!is.numeric(input$HeightPixels), 562, 
                                        ifelse(input$HeightPixels > 10000, 
                                               10000, ifelse(
                                                 input$HeightPixels < 100, 100, 
                                                 round(input$HeightPixels))
                                        )
                         ),
                         min = 100,
                         max = 10000,
                         step = 10
      )
      
      updateNumericInput(session = session, 
                         inputId = "WidthPixels", 
                         # label = tags$div(HTML(paste0('<i class="fa fa-ruler-', 
                         #                              'horizontal"></i> Plot ', 
                         #                              'Width (pixels)')
                         #                       )
                         #                  ), 
                         value = ifelse(!is.numeric(input$WidthPixels), 880, 
                                        ifelse(input$WidthPixels > 10000, 
                                               10000, ifelse(
                                                 input$WidthPixels < 100, 100, 
                                                 round(input$WidthPixels))
                                        )
                         ),
                         min = 100,
                         max = 10000,
                         step = 10
      )
      
      updateNumericInput(session = session, 
                         inputId = "log2FCFilt",
                         value = ifelse(input$log2FCFilt > 20, 20, 
                                        ifelse(input$log2FCFilt < 0, 
                                               -(input$log2FCFilt), 
                                               ifelse(
                                                 !is.numeric(input$log2FCFilt) |
                                                   is.na(input$log2FCFilt) | 
                                                   is.null(input$log2FCFilt) | 
                                                   input$log2FCFilt == "" | 
                                                   input$log2FCFilt == 0, 
                                                 2, input$log2FCFilt)
                                        )
                         ),
                         min = 0, max = 20, step = 0.1
      )
      updateNumericInput(session = session, 
                         inputId = "AdjpFilt",
                         value = ifelse(input$AdjpFilt > 1, 1, 
                                        ifelse(input$AdjpFilt < 0, 
                                               -(input$AdjpFilt), 
                                               ifelse(
                                                 !is.numeric(input$AdjpFilt) |
                                                   is.na(input$AdjpFilt) | 
                                                   is.null(input$AdjpFilt) | 
                                                   input$AdjpFilt == "" | 
                                                   input$AdjpFilt == 0, 
                                                 0.01, input$AdjpFilt)
                                        )
                         ),
                         min = 0, max = 1, step = 0.01
      )
      #print(input$AdjpFilt)
      
    })
    
    #update available samples in contrast group 1 selection buttons
    {
      #adding entries to selection based on group
      observeEvent(input$Group1Add, {
        group1_selection <- unlist(str_split(string = input$Group1Selector, 
                                             pattern = " "))
        
        #print(results_names)
        #print(group1_selection)
        if(group1_selection[2] == "Tregs"){
          available <- results_names[grepl(pattern = "cellTypeT", 
                                           x = results_names, 
                                           ignore.case = T)]
        }
        else{
          available <- results_names[grepl(pattern = "cellTypeNT", 
                                           x = results_names, 
                                           ignore.case = T)]
        }
        #print(available)
        
        if(group1_selection[1] == "Blood"){
          updateSelectInput(session = session,
                            inputId = "Group1",
                            label = "Contrast Group 1",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)]
                            ),
                            selected = c(input$Group1, available[
                              grepl(pattern = "Blood",
                                    x = available,
                                    ignore.case = T) & !(available %in%
                                                           input$Group2)])
          )
        }
        else if(group1_selection[1] == "Lymphoid"){
          updateSelectInput(session = session,
                            inputId = "Group1",
                            label = "Contrast Group 1",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)]
                            ),
                            selected = c(input$Group1, available[
                              grepl(pattern = "(LymphNode|Spleen|PeyersPatches)",
                                    x = available,
                                    ignore.case = T) & !(available %in%
                                                           input$Group2)])
          )
        }
        else if(group1_selection[1] == "Non-lymphoid"){
          updateSelectInput(session = session,
                            inputId = "Group1",
                            label = "Contrast Group 1",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)]
                            ),
                            selected = c(input$Group1, available[
                              grepl(pattern = "(Kidney|Liver|Lung|Pancreas)",
                                    x = available,
                                    ignore.case = T) & !(available %in%
                                                           input$Group2)])
          )
        }
        else if(group1_selection[1] == "Gut-associated"){
          updateSelectInput(session = session,
                            inputId = "Group1",
                            label = "Contrast Group 1",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)]
                            ),
                            selected = c(input$Group1, available[
                              grepl(pattern = "(IntraEpithelialLayer|LaminaPropriaLymphocyte)",
                                    x = available,
                                    ignore.case = T) & !(available %in%
                                                           input$Group2)])
          )
        }
        else{
          updateSelectInput(session = session,
                            inputId = "Group1",
                            label = "Contrast Group 1",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group2)]
                            ),
                            selected = c(input$Group1, available[
                              !(available %in% input$Group2)])
          )
        }
        
      })
      
      #clear selection
      observeEvent(input$Group1Clear, {
        
        updateSelectInput(session = session,
                          inputId = "Group1",
                          label = "Contrast Group 1",
                          choices = list(
                            CD4s = results_names[
                              grepl(pattern = "cellTypeNT",
                                    x = results_names,
                                    ignore.case = T) & !(results_names %in%
                                                           input$Group2)],
                            Tregs = results_names[
                              grepl(pattern = "cellTypeT",
                                    x = results_names,
                                    ignore.case = T) & !(results_names %in%
                                                           input$Group2)]
                          ),
                          selected = NULL
        )
        
      })
      
      #old button code
      {
        # #all CD4s
        # observeEvent(input$Group1SelectTconv, {
        #   
        #   updateSelectInput(session = session,
        #                     inputId = "Group1",
        #                     label = "Contrast Group 1",
        #                     choices = list(
        #                       CD4s = results_names[
        #                         grepl(pattern = "cellTypeNT",
        #                               x = results_names,
        #                               ignore.case = T) & !(results_names %in%
        #                                                      input$Group2)],
        #                       Tregs = results_names[
        #                         grepl(pattern = "cellTypeT",
        #                               x = results_names,
        #                               ignore.case = T) & !(results_names %in%
        #                                                      input$Group2)]
        #                     ),
        #                     selected = c(input$Group1, results_names[
        #                       grepl(pattern = "cellTypeNT",
        #                             x = results_names,
        #                             ignore.case = T) & !(results_names %in%
        #                                                    input$Group2)])
        #   )
        #   
        # })
        # 
        # #all Tregs
        # observeEvent(input$Group1SelectTregs, {
        #   
        #   updateSelectInput(session = session,
        #                     inputId = "Group1",
        #                     label = "Contrast Group 1",
        #                     choices = list(
        #                       CD4s = results_names[
        #                         grepl(pattern = "cellTypeNT",
        #                               x = results_names,
        #                               ignore.case = T) & !(results_names %in%
        #                                                      input$Group2)],
        #                       Tregs = results_names[
        #                         grepl(pattern = "cellTypeT",
        #                               x = results_names,
        #                               ignore.case = T) & !(results_names %in%
        #                                                      input$Group2)]
        #                     ),
        #                     selected = c(input$Group1, results_names[
        #                       grepl(pattern = "cellTypeT",
        #                             x = results_names,
        #                             ignore.case = T) & !(results_names %in%
        #                                                    input$Group2)])
        #   )
        #   
        # })
        # 
        # #all remaining
        # # observeEvent(input$Group1SelectAll, {
        # #   
        # #   updateSelectInput(session = session,
        # #                     inputId = "Group1",
        # #                     label = "Contrast Group 1",
        # #                     choices = list(
        # #                       CD4s = results_names[
        # #                         grepl(pattern = "cellTypeNT",
        # #                               x = results_names,
        # #                               ignore.case = T) & !(results_names %in%
        # #                                                      input$Group2)],
        # #                       Tregs = results_names[
        # #                         grepl(pattern = "cellTypeT",
        # #                               x = results_names,
        # #                               ignore.case = T) & !(results_names %in%
        # #                                                      input$Group2)]
        # #                     ),
        # #                     selected = c(input$Group1, results_names[
        # #                       grepl(pattern = "cellTypeNT",
        # #                             x = results_names,
        # #                             ignore.case = T) & !(results_names %in%
        # #                                                    input$Group2)], 
        # #                       results_names[
        # #                         grepl(pattern = "cellTypeT",
        # #                               x = results_names,
        # #                               ignore.case = T) & !(results_names %in%
        # #                                                      input$Group2)])
        # #   )
        # #   
        # # })
        # 
      }
      
    }
    
    #update available samples in contrast group 2 selection buttons
    {
      #adding entries to selection based on group
      observeEvent(input$Group2Add, {
        group2_selection <- unlist(str_split(string = input$Group2Selector, 
                                             pattern = " "))
        #print(group1_selection)
        if(group2_selection[2] == "Tregs"){
          available <- results_names[grepl(pattern = "cellTypeT", 
                                           x = results_names, 
                                           ignore.case = T)]
        }
        else{
          available <- results_names[grepl(pattern = "cellTypeNT", 
                                           x = results_names, 
                                           ignore.case = T)]
        }
        #print(available)
        
        if(group2_selection[1] == "Blood"){
          updateSelectInput(session = session,
                            inputId = "Group2",
                            label = "Contrast Group 2",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)]
                            ),
                            selected = c(input$Group2, available[
                              grepl(pattern = "Blood",
                                    x = available,
                                    ignore.case = T) & !(available %in%
                                                           input$Group1)])
          )
        }
        else if(group2_selection[1] == "Lymphoid"){
          updateSelectInput(session = session,
                            inputId = "Group2",
                            label = "Contrast Group 2",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)]
                            ),
                            selected = c(input$Group2, available[
                              grepl(pattern = "(LymphNode|Spleen|PeyersPatches)",
                                    x = available,
                                    ignore.case = T) & !(available %in%
                                                           input$Group1)])
          )
        }
        else if(group2_selection[1] == "Non-lymphoid"){
          updateSelectInput(session = session,
                            inputId = "Group2",
                            label = "Contrast Group 2",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)]
                            ),
                            selected = c(input$Group2, available[
                              grepl(pattern = "(Kidney|Liver|Lung|Pancreas)",
                                    x = available,
                                    ignore.case = T) & !(available %in%
                                                           input$Group1)])
          )
        }
        else if(group2_selection[1] == "Gut-associated"){
          updateSelectInput(session = session,
                            inputId = "Group2",
                            label = "Contrast Group 2",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)]
                            ),
                            selected = c(input$Group2, available[
                              grepl(pattern = "(IntraEpithelialLayer|LaminaPropriaLymphocyte)",
                                    x = available,
                                    ignore.case = T) & !(available %in%
                                                           input$Group1)])
          )
        }
        else{
          updateSelectInput(session = session,
                            inputId = "Group2",
                            label = "Contrast Group 2",
                            choices = list(
                              Tregs = results_names[
                                grepl(pattern = "cellTypeT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)],
                              Tconv = results_names[
                                grepl(pattern = "cellTypeNT",
                                      x = results_names,
                                      ignore.case = T) & !(results_names %in%
                                                             input$Group1)]
                            ),
                            selected = c(input$Group2, available[
                              !(available %in% input$Group1)])
          )
        }
        
      })
      
      #clear selection
      observeEvent(input$Group2Clear, {
        
        updateSelectInput(session = session,
                          inputId = "Group2",
                          label = "Contrast Group 2",
                          choices = list(
                            CD4s = results_names[
                              grepl(pattern = "cellTypeNT",
                                    x = results_names,
                                    ignore.case = T) & !(results_names %in%
                                                           input$Group1)],
                            Tregs = results_names[
                              grepl(pattern = "cellTypeT",
                                    x = results_names,
                                    ignore.case = T) & !(results_names %in%
                                                           input$Group1)]
                          ),
                          selected = NULL
        )
        
      })
      
      #old button code
      {
        # #all CD4s
        # observeEvent(input$Group2SelectTconv, {
        #   
        #   updateSelectInput(session = session,
        #                     inputId = "Group2",
        #                     label = "Contrast Group 2",
        #                     choices = list(
        #                       CD4s = results_names[
        #                         grepl(pattern = "cellTypeNT",
        #                               x = results_names,
        #                               ignore.case = T) & !(results_names %in%
        #                                                      input$Group1)],
        #                       Tregs = results_names[
        #                         grepl(pattern = "cellTypeT",
        #                               x = results_names,
        #                               ignore.case = T) & !(results_names %in%
        #                                                      input$Group1)]
        #                     ),
        #                     selected = c(input$Group2, results_names[
        #                       grepl(pattern = "cellTypeNT",
        #                             x = results_names,
        #                             ignore.case = T) & !(results_names %in%
        #                                                    input$Group1)])
        #   )
        #   
        # })
        # 
        # #all Tregs
        # observeEvent(input$Group2SelectTregs, {
        #   
        #   updateSelectInput(session = session,
        #                     inputId = "Group2",
        #                     label = "Contrast Group 2",
        #                     choices = list(
        #                       CD4s = results_names[
        #                         grepl(pattern = "cellTypeNT",
        #                               x = results_names,
        #                               ignore.case = T) & !(results_names %in%
        #                                                      input$Group1)],
        #                       Tregs = results_names[
        #                         grepl(pattern = "cellTypeT",
        #                               x = results_names,
        #                               ignore.case = T) & !(results_names %in%
        #                                                      input$Group1)]
        #                     ),
        #                     selected = c(input$Group2, results_names[
        #                       grepl(pattern = "cellTypeT",
        #                             x = results_names,
        #                             ignore.case = T) & !(results_names %in%
        #                                                    input$Group1)])
        #   )
        #   
        # })
        # 
        # #all remaining
        # # observeEvent(input$Group2SelectAll, {
        # #   
        # #   updateSelectInput(session = session,
        # #                     inputId = "Group2",
        # #                     label = "Contrast Group 2",
        # #                     choices = list(
        # #                       CD4s = results_names[
        # #                         grepl(pattern = "cellTypeNT",
        # #                               x = results_names,
        # #                               ignore.case = T) & !(results_names %in%
        # #                                                      input$Group1)],
        # #                       Tregs = results_names[
        # #                         grepl(pattern = "cellTypeT",
        # #                               x = results_names,
        # #                               ignore.case = T) & !(results_names %in%
        # #                                                      input$Group1)]
        # #                     ),
        # #                     selected = c(input$Group2, results_names[
        # #                       grepl(pattern = "cellTypeNT",
        # #                             x = results_names,
        # #                             ignore.case = T) & !(results_names %in%
        # #                                                    input$Group1)], 
        # #                       results_names[
        # #                         grepl(pattern = "cellTypeT",
        # #                               x = results_names,
        # #                               ignore.case = T) & !(results_names %in%
        # #                                                      input$Group1)])
        # #   )
        # #   
        # # })
        # 
      }
      
    }
    
    #enabling / disabling UI elements
    observe({
      shinyjs::toggleState(id = "GenContrasts", 
                           condition = !(length(input$Group1) < 1 | 
                                           length(input$Group2) < 1) && 
                             (is.numeric(input$log2FCFilt) &
                                is.numeric(input$AdjpFilt))
      )
    })
    #output$show_panels <- reactive({FALSE})
    #outputOptions(x = output, name = "show_panels", suspendWhenHidden = FALSE)
    shinyjs::hide(id = "res_table")
    #shinyjs::disable(id = "res_table")
    
    
    #generating contrasts via the DESeq2 results() function
    ## Generate contrasts event ----
    observeEvent(input$GenContrasts, {
      
      cat("running generate contrasts...")
      
      #splash alert for the user informing them that the analysis may take some 
      # time
      # if(!global_data$first_run_done){
      # sendSweetAlert(session = session, 
      #                title = "Computing contrasts", 
      # text = paste0("The computation may take some time during ",
      #               "which the tool may appear to not be doing ",
      #               "anything. Once the results are computed ",
      #               "they will populate the tab."),
      #                type = "info")
                  
      shinyalert(title = "Computing contrasts,\nplease wait...",
                 text = paste0("The contrasts may take some time to compute.\n",
                               "The results will populate the tab once ready."),
                 type = "info",
                 animation = F, #"slide-from-bottom",
                 inputId = "AlertFirstRunDone",
                 showCancelButton = FALSE,
                 showConfirmButton = FALSE,
                 closeOnEsc = FALSE,
                 closeOnClickOutside = FALSE)
      
      
      
      
      #   global_data$first_run_done <- T
      # }
      
      
      #ensuring that the minimum log2FC and minimum adj.p.value fields weren't 
      # left blank for the fraction of the second it takes to disable the 
      # generate contrast button
      req(input$log2FCFilt, input$AdjpFilt, input$TPCFilt)
      
      #disabling the GenContrasts button so that users may not spam it
      #is renabled at the end of this observeEvent block
      shinyjs::disable(id = "GenContrasts")
      #enabling the results table rendered object
      #shinyjs::show(id = "res_table")
      
      #resetting the previous output if it exists
      output$res_error <- NULL
      output$res_table <- NULL
      output$res_volplot <- NULL
      output$res_colplot <- NULL
      
      #print(input$TPCFilt)
      #checking which filter list to send to the contrast
      if(input$TPCFilt == "No filter"){
        tpc_filtered_list <- NULL
      }
      else{
        tpc_filtered_list <- filtered_genes_lists[[paste0("Filt", 
                                                          input$TPCFilt)]]
      }
      
      #performing the contrast
      result <- GenerateContrasts(DESeq_object = DESeq_object,
                                  group1 = input$Group1,
                                  group2 = input$Group2, 
                                  meta_data = MGI_info, 
                                  mean_norm_counts = mean_norm_counts, 
                                  tpc_filtered_list = tpc_filtered_list
      )
      #print(class(result[[2]]))
      #global_data$result <- result
      
      cat("got the result...")
      print(head(result))
      
      
      #filtering results for only genes for which actual logFC exists
      # actually bad idea because we would lose genes which are not being 
      # in one contrast group (e.g., those which are only expressed in Treg and
      # not Tconv)
      # actual_logFC <- 
      #   
      #   (is.finite(log2(rowMeans(norm_counts[,2:4])) - log2(rowMeans(norm_counts[,6:8]))))[order()]
      
      #interpretting the type of output and updating relevant UI elements
      if(length(result) == 2){
        
        #outputting error
        output$res_error <- renderText({ result[[1]] })
        
        #saving the results in a shared variable
        global_data$result <- result[[2]]
        
        result <- result[[2]]
        
      }
      else if(class(result) == "data.frame"){
        #saving the results in a shared variable
        global_data$result <- result
        
      }
      else if(class(result) == "character"){
        output$res_error <- renderText({ result })
      }
      
      #print(summary(isolate(global_data$result)))
      #calculating the fold change values for the left and right flanking 
      # columns of the NA/NaN/Inf/-Inf genes in the volcano plot
      # --can be sent to the plotting function too, but I already have that up
      # -- and running; this is for the table/data point selection and traces
      global_data$left_log <- -(ceiling(max(abs(
        isolate(global_data$result$DESeq2.Beta.Coefficients[
          is.finite(isolate(global_data$result$log2.Fold.Change))])))) + 5) #* -1
      #)))
      #print(isolate(global_data$left_log))
      global_data$right_log <- (ceiling(max(abs(
        isolate(global_data$result$DESeq2.Beta.Coefficients[
          is.finite(isolate(global_data$result$log2.Fold.Change))])))) + 5)
      #)))
      # print(isolate(global_data$right_log))
      
      
      #processing data frame to be displayed as DataTable
      if(class(result) == "data.frame"){
        
        #replacing zero p.values/adj.p.values with smallest IEEE supported value
        global_data$result$p.value[global_data$result$p.value == 0] <- 
          .Machine$double.xmin
        global_data$result$Adjusted.p.value[
          global_data$result$Adjusted.p.value == 0] <- .Machine$double.xmin
        
        result <- global_data$result
        
        #removing extra columns
        result$Standard.Error <- NULL
        result$Statistic <- NULL
        result$p.value <- NULL
        
        #formatting the numeric columns out here rather than through DT::format*
        # so that the filter ranges do not end up with 15 decimal places
        result$Mean.Counts <- round(x = result$Mean.Counts, digits = 3)
        result$log2.Fold.Change <- round(x = result$log2.Fold.Change,
                                         digits = 3)
        result$DESeq2.Beta.Coefficients <- 
          round(x = result$DESeq2.Beta.Coefficients, digits = 3)
        result$Adjusted.p.value <- signif(x = result$Adjusted.p.value,
                                          digits = 3)
        
        #adding selection column for sorting selected rows -- not being used atm
        result$Selected <- "No"
        
        #outputting results
        output$res_table <- DT::renderDataTable({ 
          DT::datatable(elementId = "res_table", 
                        data = result, 
                        #caption = "Differential expression results",
                        colnames = gsub(pattern = "\\.", 
                                        replacement = " ", 
                                        x = colnames(result), 
                                        ignore.case = T),
                        #extensions = "Buttons",
                        options = list(searching = T, 
                                       searchHighlight = TRUE,
                                       dom = "lrtip",
                                       #buttons = c('csv', 'excel', 'pdf'),
                                       #order = list(list(5, 'asc'))
                                       order = list(list(6, 'asc'))
                        ), 
                        filter = list(position = "top", 
                                      plain = F
                        ),
                        editable = F, 
                        rownames = F) #%>% 
          # DT::formatRound(columns = c("Mean.Counts"),
          #                 digits = 2,
          #                 mark = "") %>%
          # DT::formatRound(columns = c("log2.Fold.Change"),
          #                 digits = 3,
          #                 mark = "") %>%
          # DT::formatSignif(columns = c("Adjusted.p.value"),
          #                  digits = 3,
          #                  mark = "")
        })
        
        #used for updating DataTable values for the 'Selected' column
        global_data$dt <- result
        
      }
      
      cat("got all the way down here to line 2274...")
      
      
      #checking if a contrast was generated to create the volcano plot
      if(class(result) == "data.frame"){
        
        #interpreting the label names for the plot title
        names_group1 <- names(results_names[results_names %in% 
                                              isolate(input$Group1)])
        names_group2 <- names(results_names[results_names %in% 
                                              isolate(input$Group2)])
        
        #turning off warning to hide the persistent "Warning: 'scattergl' 
        # objects don't have these attributes: 'hoveron'" warnings. Comment the
        # below code when debugging the plotting code
        #restoreWarn <- getOption("warn")
        #options(warn = -1) 
        
        #generating the ggplotly plot
        # output$res_volplot <- renderPlotly({
        #   GenerateVolcanoPlot(plot_data = isolate(global_data$result),
        #                       names_group1 = names_group1,
        #                       names_group2 = names_group2, 
        #                       format = tolower(isolate(input$SavePlotFormat)), 
        #                       width = isolate(input$WidthPixels), 
        #                       height = isolate(input$HeightPixels)
        #                       )
        p <- GenerateVolcanoPlot(plot_data = isolate(global_data$result),
                                 names_group1 = names_group1,
                                 names_group2 = names_group2, 
                                 format = tolower(isolate(input$SavePlotFormat)), 
                                 width = isolate(input$WidthPixels), 
                                 height = isolate(input$HeightPixels), 
                                 log2FCFilt = input$log2FCFilt, 
                                 AdjpFilt = input$AdjpFilt
        )
        output$res_volplot <- renderPlotly({p})
        p <- plotly_build(p)
        global_data$base_traces <- length(p$x$data)
        global_data$trace_count <- global_data$base_traces-1
        #plotly traces are zero indexed so trace_count is base_traces-1
        
        #restoring warnings
        #options(warn = restoreWarn) 
        
      }
      
      #renabling the disabled GenContrasts button
      shinyjs::enable(id = "GenContrasts")
      #shinyjs::enable(id = "res_table")
      
      #showing hidden ui controls and elements
      shinyjs::show(id = "res_table")
      output$show_panels <- reactive({TRUE}) #for multiple conditionalPanels
      outputOptions(x = output, name = "show_panels", suspendWhenHidden = FALSE)
      
      # #update the SavePlotFormat once more to allow saving the default SVG
      # updateRadioGroupButtons(session = session, 
      #                         inputId = "SavePlotFormat", 
      #                         selected = isolate(input$SavePlotFormat))
      
      #resuming the DT row selection observer
      #res_table_row_observer$resume()
      
      #remove popup shinyalert from earlier using second alert with immediate 
      # set to TRUE
      closeAlert()
      # shinyalert(title = "Computation Complete", 
      #            text = paste0("Results are being populated"), 
      #            type = "info", 
      #            animation = F, #"slide-from-bottom", 
      #            inputId = "RemoveAlert", 
      #            immediate = TRUE, timer = 1, 
      #            showCancelButton = FALSE, 
      #            showConfirmButton = FALSE, 
      #            closeOnEsc = FALSE, 
      #            closeOnClickOutside = FALSE)
      
    })
    ## end of generate contrasts observeEvent ----
    
    #generating bar plot of normalised counts for the genes of the rows selected
    # in the rendered DataTable
    observeEvent(input$res_table_rows_selected, ignoreNULL = F, {
      
      #resetting any existing plot
      output$res_colplot <- NULL
      
      #remove any existing plotly traces on the volcano plot
      # plotly traces are zero indexed so min count is global_data$base_traces-1
      #if(global_data$trace_count > 2){
      if(global_data$trace_count > global_data$base_traces-1){
        plotlyProxy(outputId = "res_volplot", session = session) %>%
          plotlyProxyInvoke("deleteTraces",
                            #as.list(c(2:global_data$trace_count))
                            as.list(c(
                              global_data$base_traces:global_data$trace_count))
          )
        
        global_data$trace_count <- global_data$base_traces-1
        #print(global_data$trace_count)
      }
      
      #if there is at least one row selected
      if(!is.null(input$res_table_rows_selected)){
        
        #hard limiting max selectable rows to first 50 selected
        if(length(input$res_table_rows_selected) > 50){
          selected_rows <- input$res_table_rows_selected[1:50]
          DT::dataTableProxy("res_table") %>% selectRows(selected_rows)
        }
        else{
          selected_rows <- input$res_table_rows_selected
        }
        
        selection <- isolate(global_data$result[selected_rows, 1:2])
        
        #generating the counts plot
        output$res_colplot <- renderPlotly({
          GenerateColPlot(mean_norm_counts = mean_norm_counts, 
                          selection = selection,
                          format = tolower(isolate(input$SavePlotFormat)), 
                          width = isolate(input$WidthPixels), 
                          height = isolate(input$HeightPixels)
          )
        })
        
        #adding traces to the plotly volcano plot
        #need cols 5 and 9 from selection
        selection <- isolate(global_data$result[selected_rows,])
        #print(summary(selection))
        
        #colnames: 1-Gene.ID 2-Gene.Symbol 3-Description 4-Mean.Counts 
        #           5-log2.Fold.Change 6-DESeq2.Beta.Coefficients 7-p.value 
        #           8-Adjusted.p.value 9-Standard.Error 10-Statistic
        
        #output$debug_console <- renderPrint(selection)
        
        #plotly does not show traces for single items -- 
        #duplicating item if single
        if(nrow(selection) == 1){
          #checking if log2.Fold.Change exists, then use DESeq2.Beta.Coefficients
          # else use max/min for the right/left columns
          if(is.finite(selection[1,5])){
            selected_x <- c(selection[1,6],selection[,6])
          }
          else{
            if(selection[1,6] < 0) selected_x <- 
                c(isolate(global_data$left_log), isolate(global_data$left_log))
            else if(selection[1,6] > 0) selected_x <- 
                c(isolate(global_data$right_log), isolate(global_data$right_log))
          }
          #selected_y <- c(selection[1,9],selection[,9])
          selected_y <- c(selection[1,8],selection[,8])
          
          #if MGI or other symbols are missing
          #if(selection[,2] == ""){
          if(is.na(selection[,2]) | selection[,2] == ""){
            # selected_text <- c(str_extract(selection[,1], '\\w*')[1],
            #                    str_extract(selection[,1], '\\w*'))
            selected_text <- c(gsub(pattern = ";.*", 
                                    replacement = "", 
                                    x = selection[,1], 
                                    ignore.case = T)[1],
                               gsub(pattern = ";.*", 
                                    replacement = "", 
                                    x = selection[,1], 
                                    ignore.case = T)
            )
            #print(selected_text)
          }
          else{
            # selected_text <- c(str_extract(selection[,2], '\\w*')[1],
            #                    str_extract(selection[,2], '\\w*'))
            selected_text <- c(gsub(pattern = ";.*", 
                                    replacement = "", 
                                    x = selection[,2], 
                                    ignore.case = T)[1],
                               gsub(pattern = ";.*", 
                                    replacement = "", 
                                    x = selection[,2], 
                                    ignore.case = T)
            )
            #print(selected_text)
          }
        }
        else{
          #selected_x <- selection[,5]
          selected_x <- selection[,6]
          #checking if log2.Fold.Change exists, then use DESeq2.Beta.Coefficients
          # else use max/min for the right/left columns
          selected_x[!is.finite(selection[,5]) & selection[,6] < 0] <- 
            isolate(global_data$left_log)
          selected_x[!is.finite(selection[,5]) & selection[,6] > 0] <- 
            isolate(global_data$right_log)
          
          #selected_y <- selection[,9]
          selected_y <- selection[,8]
          
          #if MGI or other symbols are missing
          #if(selection[,2] == ""){
          if(any(is.na(selection[,2])) | any(selection[,2] == "")){
            #print(selection[,2])
            #selected_text <- str_extract(selection[,1], '\\w*')
            selected_text <- gsub(pattern = ";.*", 
                                  replacement = "", 
                                  x = selection[,1], 
                                  ignore.case = T)
            selected_text[-(which(is.na(selection[,2]) | 
                                    selection[,2] == ""))] <- 
              gsub(pattern = ";.*", 
                   replacement = "", 
                   x = selection[-(which(is.na(selection[,2]) | 
                                           selection[,2] == "")),2], 
                   ignore.case = T)
            # selected_text[-(which(is.na(selection[,2]) | 
            #                         selection[,2] == ""))] <- 
            #   str_extract(selection[-(which(is.na(selection[,2]) | 
            #                                   selection[,2] == "")),2], '\\w*')
            #print(selected_text)
          }
          else{
            #selected_text <- str_extract(selection[,2], '\\w*')
            selected_text <- gsub(pattern = ";.*", 
                                  replacement = "", 
                                  x = selection[,2], 
                                  ignore.case = T)
          }
        }
        
        plotlyProxy(outputId = "res_volplot", session = session) %>%
          plotlyProxyInvoke("addTrace",
                            list(#x = c(selection[1,5],selection[,5]),
                              #y = -log10(c(selection[1,9], selection[,9])),
                              x = selected_x,
                              y = -log10(selected_y),
                              type = "scattergl",
                              name = "Selected genes",
                              mode = "markers+text",
                              marker = list(#symbol = "x",
                                symbol = "circle-open", 
                                size = 10, 
                                #colorscale = "Viridis"
                                color = "blue"
                              ), 
                              # text = c(str_extract(selection[,2], '\\w*')[1],
                              #          str_extract(selection[,2], '\\w*')), 
                              text = selected_text,
                              textposition = "top right" ,
                              #textfont_size = 4, 
                              textfont = list(family = "sans serif", 
                                              size = 10, 
                                              color = "blue"
                              )
                            )
          )
        
        # p <- plotlyProxy(outputId = "res_volplot", session = session)$x$data
        # saveRDS(object = p, file = "p_proxy_build.rds")
        
        global_data$trace_count <- global_data$trace_count + 1
        #print(global_data$trace_count)
        #output$debug_console <- renderPrint(global_data$trace_count)
        
        #output$debug_console <- renderPrint(selection)
        
        #update the selection column of the selected rows
        #was designed so that selected rows can then be sorted on to move them to
        #the top. however, the code currently invokes a perpetual looping update 
        #of the table and doesn't allow the selection of more than one row either
        #probably due to reactivity of the row selection event -- maybe try 
        #binding to an action button instead of automatic updates via row 
        #selection? ---- code is now fixed
        changes <- data.frame(row = selected_rows,
                              #col = 6, #column index of "selected" (index starts at 0)
                              col = 7, 
                              value = "Yes")
        
        DT::dataTableProxy("res_table") %>% 
          DT::editData(data = isolate(global_data$dt), 
                       info = changes,
                       rownames = F, 
                       resetPaging = F, 
                       clearSelection = "none")
        
      }
      #else if the last row has been deselected, then reset selection column
      else if(!is.null(input$res_table_row_last_clicked)){
        
        changes <- data.frame(row = isolate(input$res_table_row_last_clicked),
                              #col = 6, #column index of "selected" (index starts at 0)
                              col = 7, 
                              value = "No")
        
        DT::dataTableProxy("res_table") %>% 
          DT::editData(data = isolate(global_data$dt), 
                       info = changes,
                       rownames = F, 
                       resetPaging = F, 
                       clearSelection = "all")
      }
      
    })
    
    
    #resetting rows selected in the results table
    observeEvent(input$ResetTable, {
      proxy <- dataTableProxy("res_table")
      
      proxy %>% selectRows(NULL)
      
      output$res_colplot <- NULL
      
      #removing traces from the plotly volcano plot
      # plotly traces are zero indexed so min count is global_data$base_traces-1
      if(global_data$trace_count > global_data$base_traces-1){
        plotlyProxy(outputId = "res_volplot", session = session) %>%
          plotlyProxyInvoke("deleteTraces",
                            #as.list(c(2:global_data$trace_count))
                            as.list(c(
                              global_data$base_traces:global_data$trace_count))
          )
        
        global_data$trace_count <- global_data$base_traces-1
        #print(global_data$trace_count)
      }
      # if(global_data$trace_count > 2){
      #   plotlyProxy(outputId = "res_volplot", session = session) %>% 
      #     plotlyProxyInvoke("deleteTraces", 
      #                       #as.list(c(2:global_data$trace_count))
      #                       as.list(c(3:global_data$trace_count))
      #                       # as.list(c(
      #                       #   (global_data$trace_count-1):global_data$trace_count))
      #     )
      #   
      #   global_data$trace_count <- 2
      # }
      
    })
    
    
    #download the complete results
    output$DownloadTableAll <- downloadHandler(filename = "Results.tab", 
                                               content = function(file){
                                                 write.table(
                                                   x = global_data$result, 
                                                   file = file, 
                                                   quote = F, 
                                                   sep = "\t", 
                                                   col.names = T, 
                                                   row.names = F)
                                               })
    
    #download the filtered results
    output$DownloadTableFilt <- downloadHandler(
      filename = "Results_filtered.tab", 
      content = function(file){
        write.table(x = 
                      global_data$result[input$res_table_rows_all, , drop = F], 
                    file = file, 
                    quote = F, 
                    sep = "\t", 
                    col.names = T, 
                    row.names = F)
      })
    
    #download the selected rows only 
    output$DownloadTableSelected <- downloadHandler(
      filename = "Results_selected.tab", 
      content = function(file){
        write.table(x = global_data$result[
          input$res_table_rows_selected, , drop = F], 
          file = file, 
          quote = F, 
          sep = "\t", 
          col.names = T, 
          row.names = F)
      })
    
    #(de)select rows in the results table by (de)selecting it in volcano plot
    observe({
      
      #the req() circumvents the whole exists() and !is.null requirements in the
      #outer if-else block below, reducing the check to just the 50 hard limit
      req(event_data("plotly_click", 
                     source = "volplot", 
                     priority = "event"))
      
      volplot_selected <- event_data("plotly_click", 
                                     source = "volplot", 
                                     priority = "event")
      proxy <- dataTableProxy("res_table")
      
      already_selected <- isolate(input$res_table_rows_selected)
      # proxy %>% selectRows(c(volplot_selected$key, already_selected))
      
      #the following code is designed to allow deselection, but it is bugged for
      #some reason and requires the user to click the last selected node twice
      #to deselect it will need to work on this in later updates if required as
      #a feature ----- the bug is resolved by 'priority = "event"' above
      #Also, hard limiting to max selections to 50
      # if(exists("volplot_selected") & !is.null(volplot_selected) &
      #    length(already_selected) < 50){
      if(length(already_selected) < 50){
        
        if(length(intersect(volplot_selected$key, already_selected)) > 0){
          #print("entering 'if' block")
          if(length(already_selected) > 1){
            #print("entering 'if>if' block")
            proxy %>% selectRows(as.numeric(
              setdiff(already_selected, volplot_selected$key)))
          }
          else{
            #print("entering 'if>else' block")
            proxy %>% selectRows(NULL)
          }
        }
        else if(length(union(volplot_selected$key, already_selected)) > 0){
          #print("entering 'else if' block")
          proxy %>% selectRows(as.numeric(
            union(volplot_selected$key, already_selected)))
        }
        # else if(length(volplot_selected$key) > 0 &
        #         length(already_selected) < 1){
        #   proxy %>% selectRows(as.numeric(volplot_selected$key))
        # }
        # else{
        #   print("entering 'else' block")
        #   proxy %>% selectRows(NULL)
        # }
        
      }
      
      #volplot_selected <- NULL
      #print(isolate(global_data$trace_count))
      
    })
    
    
    #creating reactive expression to group multiple inputs for the observeEvent
    plotFormatAndSize <- reactive({
      list(input$SavePlotFormat, input$WidthPixels, input$HeightPixels)
    })
    #change the snapshot save format and size of plotly plots
    observeEvent(plotFormatAndSize(), {
      
      #formatting the input
      format <- tolower(input$SavePlotFormat)
      width <- input$WidthPixels
      height <- input$HeightPixels
      
      #updating the format for the volcano plot
      plotlyProxy(outputId = "res_volplot", session = session) %>%
        plotlyProxyInvoke("reconfig",
                          displaylogo = F, 
                          displayModeBar = T, 
                          toImageButtonOptions = 
                            list(format = format, 
                                 width = width, 
                                 height = height
                            ),
                          modeBarButtonsToRemove = list(
                            "sendDataToCloud", "zoom2d", "zoomIn2d",
                            "zoomOut2d", "pan2d", "select2d", "lasso2d",
                            "autoScale2d", "hoverClosestCartesian",
                            "hoverCompareCartesian"
                          )
        )
      
      #updating the format for the counts plot
      plotlyProxy(outputId = "res_colplot", session = session) %>%
        plotlyProxyInvoke("reconfig",
                          displaylogo = F, 
                          displayModeBar = T, 
                          toImageButtonOptions = 
                            list(format = format, 
                                 width = width, 
                                 height = height
                            ),
                          modeBarButtonsToRemove = list(
                            "sendDataToCloud", "zoom2d", "zoomIn2d",
                            "zoomOut2d", "pan2d", "select2d", "lasso2d",
                            "autoScale2d", "toggleSpikelines"#, 
                            #"hoverClosestCartesian",
                            #"hoverCompareCartesian"
                          )
        )
      
    })
    
    
    #display PCA, t-SNE and UMAP plots
    observeEvent(input$TriplePlot, {
      
      output$PCAplot1 <- NULL
      output$PCAplot2 <- NULL
      output$PCAplot3 <- NULL
      
      # output$debug_console <- renderPrint(input$TriplePlot == "PCA")
      
      if(input$TriplePlot == "PCA"){
        #output$debug_console <- renderText(class(input$TriplePlot))
        #output$PCAplot1 <- renderPlotly({
        output$PCAplot1 <- renderPlot({
          #ggplotly(PCA_plots[[1]]) %>% toWebGL()
          PCA_plots[[1]]
        })
        output$PCAplot2 <- renderPlot({
          PCA_plots[[2]]
        })
        output$PCAplot3 <- renderPlot({
          PCA_plots[[3]]
        })
      }
      else if(input$TriplePlot == "t-SNE"){
        output$PCAplot1 <- renderPlot({
          tSNE_plots[[1]]
        })
        output$PCAplot2 <- renderPlot({
          tSNE_plots[[2]]
        })
        output$PCAplot3 <- renderPlot({
          tSNE_plots[[3]]
        })
      }
      else if(input$TriplePlot == "UMAP"){
        output$PCAplot1 <- renderPlot({
          UMAP_plots[[1]]
        })
        output$PCAplot2 <- renderPlot({
          UMAP_plots[[2]]
        })
        output$PCAplot3 <- renderPlot({
          UMAP_plots[[3]]
        })
      }
      else if(input$TriplePlot == "PCA - tissue means"){
        output$PCAplot1 <- renderPlot({
          PCA_means_plots[[1]]
        })
        output$PCAplot2 <- renderPlot({
          PCA_means_plots[[2]]
        })
        output$PCAplot3 <- renderPlot({
          PCA_means_plots[[3]]
        })
      }
      
    })
    
    
    #display TCR chord diagrams
    observeEvent(input$TCRChordSelect, {
      output$TCRChordImage <- renderImage(expr = {
        
        img <- readPNG(
          source = paste0("www/TCR_plots/", 
                          TCR_plots_list[grep(pattern = input$TCRChordSelect, 
                                              x = TCR_plots_list, 
                                              ignore.case = T)]), 
          native = T, info = T)
        
        list(src = paste0("www/TCR_plots/", 
                          TCR_plots_list[grep(pattern = input$TCRChordSelect, 
                                              x = TCR_plots_list, 
                                              ignore.case = T)]),
             contentType = 'image/png',
             width = dim(img)[2]*0.15,
             height = dim(img)[1]*0.15#,
             #alt = "This is alternate text"
        )
      }, deleteFile = F)
    })
    
    
    tsne_aging_img <- function(filename){
      img <- readPNG(source = filename, native = T, info = T)
      list(
        src = filename, 
        contentType = 'image/png', 
        width = dim(img)[2]*0.4, 
        height = dim(img)[1]*0.4
      )
    }
    
    observeEvent(input$browser, browser())
    
    # aging tSNE plots ----
    observeEvent(input$agingTSNEtissue, {
      tissue_files <- paste0(
        "www/Aging_tsne_plots/", 
        aging_tsne_plots_list[
          grep(
            pattern = input$agingTSNEtissue, 
            x = aging_tsne_plots_list, 
            ignore.case = T)
          ]
      )
      # ordering them manually - they're read in in this order due to the numbers:
      # [1] "xx100-week__cluster.png"
      # [2] "xx12-week__cluster.png" 
      # [3] "xx20-week__cluster.png" 
      # [4] "xx30-week__cluster.png" 
      # [5] "xx52-week__cluster.png" 
      # [6] "xx8-week__cluster.png" 
      tissue_files <- tissue_files[c(6,2,3,4,5,1)]
      
      output$agingTSNEImage1 <- renderImage(
        expr = tsne_aging_img(tissue_files[1]), deleteFile = F
      )
      output$agingTSNEImage2 <- renderImage(
        expr = tsne_aging_img(tissue_files[2]), deleteFile = F
      )
      output$agingTSNEImage3 <- renderImage(
        expr = tsne_aging_img(tissue_files[3]), deleteFile = F
      )
      output$agingTSNEImage4 <- renderImage(
        expr = tsne_aging_img(tissue_files[4]), deleteFile = F
      )
      output$agingTSNEImage5 <- renderImage(
        expr = tsne_aging_img(tissue_files[5]), deleteFile = F
      )
      output$agingTSNEImage6 <- renderImage(
        expr = tsne_aging_img(tissue_files[6]), deleteFile = F
      )
      output$agingTSNETitle1 <- renderText(paste("t-SNE", input$agingTSNEtissue, "8 weeks"))
      output$agingTSNETitle2 <- renderText(paste("t-SNE", input$agingTSNEtissue, "12 weeks"))
      output$agingTSNETitle3 <- renderText(paste("t-SNE", input$agingTSNEtissue, "20 weeks"))
      output$agingTSNETitle4 <- renderText(paste("t-SNE", input$agingTSNEtissue, "30 weeks"))
      output$agingTSNETitle5 <- renderText(paste("t-SNE", input$agingTSNEtissue, "52 weeks"))
      output$agingTSNETitle6 <- renderText(paste("t-SNE", input$agingTSNEtissue, "100 weeks"))
    })
    
    
    tsne_microbiome_img <- function(filename){
      img <- readPNG(source = filename, native = T, info = T)
      list(
        src = filename, 
        contentType = 'image/png', 
        width = dim(img)[2]*0.4, 
        height = dim(img)[1]*0.4
      )
    }
    
    
    # microbiome tSNE plots ----
    observeEvent(input$microbiomeTSNEtissue, {
      tissue_files <- paste0(
        "www/Microbiome_tsne_plots/", 
        microbiome_tsne_plots_list[grep(pattern = input$microbiomeTSNEtissue, 
                                        x = microbiome_tsne_plots_list, 
                                        ignore.case = TRUE)])

      output$microbiomeTSNEImage1 <- renderImage(
        expr = tsne_microbiome_img(tissue_files[1]), deleteFile = F
      )
      output$microbiomeTSNEImage2 <- renderImage(
        expr = tsne_microbiome_img(tissue_files[2]), deleteFile = F
      )
      output$microbiomeTSNEImage3 <- renderImage(
        expr = tsne_microbiome_img(tissue_files[3]), deleteFile = F
      )

      output$microbiomeTSNETitle1 <- renderText(paste("t-SNE", input$microbiomeTSNEtissue, "Cohoused"))
      output$microbiomeTSNETitle2 <- renderText(paste("t-SNE", input$microbiomeTSNEtissue, "Gnotobiotic"))
      output$microbiomeTSNETitle3 <- renderText(paste("t-SNE", input$microbiomeTSNEtissue, "SPF"))

    })
    
    
    
    #display Pathview KEGG Pathways
    observeEvent(input$PathviewPlotSelect, {
      output$PathviewPlotImage <- renderImage(expr = {
        
        pathway_id <- gsub(pattern = " - .*", 
                           replacement = "", 
                           x = input$PathviewPlotSelect, 
                           ignore.case = T)
        
        img <- readPNG(
          source = paste0("www/Pathview_plots/", 
                          Pathview_plots_list[grep(pattern = pathway_id, 
                                                   x = Pathview_plots_list, 
                                                   ignore.case = T)]), 
          native = T, info = T)
        
        list(src = paste0("www/Pathview_plots/", 
                          Pathview_plots_list[grep(pattern = pathway_id, 
                                                   x = Pathview_plots_list, 
                                                   ignore.case = T)]),
             contentType = 'image/png',
             width = dim(img)[2],#*0.35,
             height = dim(img)[1]#*0.35#,
             #alt = "This is alternate text"
        )
      }, deleteFile = F)
    })
    
    #display Pathview Legend and Legend2
    output$PathviewPlotLegend <- renderImage(expr = {
      
      img <- readPNG(source = "www/Legend.png", 
                     native = T, info = T)
      
      list(src = "www/Legend.png",
           contentType = 'image/png',
           width = dim(img)[2],
           height = dim(img)[1]#,
           #alt = "This is alternate text"
      )
    }, deleteFile = F)
    output$PathviewPlotLegend2 <- renderImage(expr = {
      
      img <- readPNG(source = "www/Legend2.png", 
                     native = T, info = T)
      
      list(src = "www/Legend2.png",
           contentType = 'image/png',
           width = dim(img)[2]*0.6,
           height = dim(img)[1]*0.6#,
           #alt = "This is alternate text"
      )
    }, deleteFile = F)
    
    #download Pathview pathways
    output$DownloadPathviewSelected <- downloadHandler(
      filename = function() {
        paste0(input$PathviewPlotSelect, ".zip")
      },
      #filename = paste0(input$PathviewPlotSelect, ".zip"),
      content = function(file){
        
        pathway_id <- gsub(pattern = " - .*",
                           replacement = "",
                           x = input$PathviewPlotSelect,
                           ignore.case = T)
        
        zip(zipfile = file, 
            files = c(paste0("www/Pathview_plots/",
                             Pathview_plots_list[grep(pattern = pathway_id,
                                                      x = Pathview_plots_list,
                                                      ignore.case = T)]), 
                      "www/Legend.png", "www/Legend2.png"), 
            flags = "-j")
        
      }, 
      contentType = "application/zip")

    output$DownloadPathviewAll <- downloadHandler(
      filename = function() {
        paste0("all_pathways.zip")
      },
      #filename = "all_pathways.zip",
      content = function(file){
        
        zip(zipfile = file, 
            files = c(paste0("www/Pathview_plots/", Pathview_plots_list), 
                      "www/Legend.png", "www/Legend2.png"), 
            flags = "-j")
        
      },
      contentType = "application/zip")
    
    # download aging tSNE plots ----
    output$DownloadAgingTSNEtissue <- downloadHandler(
      
      filename = function() {
        paste0(input$agingTSNEtissue, "_aging_tSNE.zip")
      },
      content = function(file){

        zip(
          zipfile = file,
          files = 
            paste0(
              "www/Aging_tsne_plots/", 
              aging_tsne_plots_list[
                grep(pattern = input$agingTSNEtissue, x = aging_tsne_plots_list, ignore.case = T)
              ]
            ),
          flags = "-j"
        )
      }, contentType = "application/zip"
    )
    
    
    output$DownloadAgingTSNEAll <- downloadHandler(
      filename = function() {
        paste0("all_aging_tSNE.zip")
      },
      content = function(file){
        
        zip(zipfile = file, 
            files = paste0("www/Aging_tsne_plots/", aging_tsne_plots_list), 
            flags = "-j")
        
      },
      contentType = "application/zip"
    )
    
    # download microbiome tSNE plots ----
    output$DownloadMicrobiomeTSNEtissue <- downloadHandler(
      
      filename = function() {
        paste0(input$microbiomeTSNEtissue, "_microbiome_tSNE.zip")
      },
      content = function(file){
        
        zip(
          zipfile = file,
          files = 
            paste0(
              "www/Microbiome_tsne_plots/", 
              microbiome_tsne_plots_list[
                grep(pattern = input$microbiomeTSNEtissue, x = microbiome_tsne_plots_list, ignore.case = T)
              ]
            ),
          flags = "-j"
        )
      }, contentType = "application/zip"
    )
    
    
    output$DownloadMicrobiomeTSNEAll <- downloadHandler(
      filename = function() {
        paste0("all_microbiome_tSNE.zip")
      },
      content = function(file){
        
        zip(zipfile = file, 
            files = paste0("www/Microbiome_tsne_plots/", microbiome_tsne_plots_list), 
            flags = "-j")
        
      },
      contentType = "application/zip"
    )
    
    
    # download TCR plots ----
    output$DownloadTCRSelected <- downloadHandler(
      filename = function() {
        paste0(input$TCRChordSelect, "_chord_diagram.png")
      },
      #filename = paste0(input$TCRChordSelect, "_chord_diagram.png"),
      content = function(file){
        
        file.copy(paste0("www/TCR_plots/",
                         TCR_plots_list[grep(pattern = input$TCRChordSelect,
                                             x = TCR_plots_list,
                                             ignore.case = T)]),
                  file,
                  copy.mode = F,
                  copy.date = F
        )
      },
      contentType = "image/png")
    
    output$DownloadTCRAll <- downloadHandler(
      filename = function() {
        paste0("all_chord_diagrams.zip")
      },
      #filename = "all_chord_diagrams.zip",
      content = function(file){
        
        zip(zipfile = file, 
            files = paste0("www/TCR_plots/", TCR_plots_list), 
            flags = "-j")
        
      },
      contentType = "application/zip")
    
    #downloading Supplementary Materials
    output$MultiQCBulk <- downloadHandler(
      filename = function() {
        paste0("QC and pre-processing report.zip")
      },
      content = function(file){
        file.copy("www/MultiQC_report/QC and pre-processing report.zip",
                  file,
                  copy.mode = F,
                  copy.date = F
        )
      },
      contentType = "application/zip")
    output$GlobalBulk <- downloadHandler(
      filename = function() {
        paste0("PCA, t-SNE and UMAP plots.zip")
      },
      content = function(file){
        zip(zipfile = file,
            files = paste0("www/Global_bulk/", 
                           list.files(path = "www/Global_bulk/")),
            flags = "-j"
        )
      },
      contentType = "application/zip")
    #"GSEABulk", DiffExprsBulk
    output$DiffExprsBulk <- downloadHandler(
      filename = function() {
        paste0("Counts and differential expression table.xlsx")
      },
      content = function(file){
        file.copy("www/DESeq2_table/output_table.xlsx",
                  file,
                  copy.mode = F,
                  copy.date = F
        )
      },
      contentType = "xlsx")
    output$GSEABulk <- downloadHandler(
      filename = function() {
        paste0("GSEA table.xlsx")
      },
      content = function(file){
        file.copy("www/GSEA_table/output_gage.xlsx",
                  file,
                  copy.mode = F,
                  copy.date = F
        )
      },
      contentType = "xlsx")
    
    
    #dynamic tutorials for rintrojs
    {
      steps <- reactive(data.frame(element = c("#Group1 + .selectize-control", 
                                               "#Group1Selector + .selectize-control", 
                                               "#Group1Add", 
                                               "#Group1Clear", 
                                               "#Group2 + .selectize-control", 
                                               "#AddOps", 
                                               "#TPCFilt + .selectize-control", 
                                               "#log2FCFilt", 
                                               "#AdjpFilt", 
                                               "#GenContrasts", 
                                               #"#res_table", 
                                               #"#DownloadTableAll", 
                                               #"#DownloadTableFilt", 
                                               #"#DownloadTableSelected", 
                                               #"#ResetTable", 
                                               "#mydropdown"),
                                   intro = c(paste0("Start by selecting the ", 
                                                    "tissue(s) that need to be ", 
                                                    "contrasted. Selected ", 
                                                    "tissue(s) can be removed ", 
                                                    "by clicking on the name ", 
                                                    "and pressing the 'backspace' ", 
                                                    "or 'delete' key."), 
                                             paste0("Groups of tissues can be ", 
                                                    "selected in bulk from ", 
                                                    "here..."),
                                             paste0("...and added using this ", 
                                                    "button."),
                                             paste0("The whole selection can ", 
                                                    "be cleared by clicking ", 
                                                    "this button."),
                                             paste0("Tissue(s) to be contrasted ", 
                                                    "against can be selected ", 
                                                    "here. Note that tissues ", 
                                                    "already selected in the ", 
                                                    "first group are not ", 
                                                    "available here and vice ", 
                                                    "versa."), 
                                             paste0("This checkbox shows and ", 
                                                    "hides additional options ", 
                                                    "for the analysis."),
                                             paste0("Transcripts per cell ", 
                                                    "filtering filters out low ", 
                                                    "expression genes from the ", 
                                                    "results and plots. The ", 
                                                    "filtering criteria is at ", 
                                                    "least the selected amount ", 
                                                    "of counts per cell in at ", 
                                                    "least one Treg tissue."), 
                                             paste0("The minimum absolute log2 ", 
                                                    "fold change for a gene to ", 
                                                    "be considered ", 
                                                    "differentially expressed."),
                                             paste0("The maximum adjusted p.", 
                                                    "value for a gene to be ", 
                                                    "considered differentially ", 
                                                    "expressed."),
                                             paste0("Clicking this button will ", 
                                                    "generate the contrasts ", 
                                                    "between the groups ", 
                                                    "selected above. This ", 
                                                    "button remains disabled ", 
                                                    "till both groups have at ", 
                                                    "least one tissue selected."),
                                             #paste0(),
                                             paste0("This button can be ", 
                                                    "clicked to show a summary ", 
                                                    "of the instructions ", 
                                                    "and features."))))
      
      observeEvent(input$TutorialButton,{
        updateCheckboxInput(session = session, 
                            inputId = "AddOps", 
                            value = TRUE)
        
        introjs(session = session , options = list(steps=steps()))
      })
      
      steps2 <- reactive(data.frame(element = c("#res_table", 
                                                "#DownloadTableAll", 
                                                "#DownloadTableFilt", 
                                                "#DownloadTableSelected", 
                                                "#ResetTable", 
                                                "#PlotSaveOps", 
                                                "#SavePlotFormat", 
                                                "#WidthPixels", 
                                                "#HeightPixels", 
                                                "#res_volplot" 
      ),
      intro = c(paste0("The results of the ", 
                       "contrast are tabulated ", 
                       "here. Column headers can ", 
                       "be clicked to sort by ", 
                       "that column (sorting can ", 
                       "be done on multiple ", 
                       "columns by pressing the ", 
                       "'shift' key and clicking ", 
                       "on the other header). ", 
                       "Columns can be searched ", 
                       "or filtered by clicking ", 
                       "on the blank fields ", 
                       "under the column headers.", 
                       " Clicking a row will ", 
                       "select the row in the ", 
                       "table and highlight the ", 
                       "respective gene in the ", 
                       "volcano plot (upto 50 ", 
                       "rows can be selected by ", 
                       "pressing the 'shift' or ", 
                       "'control' key)."), 
                paste0("Download the complete ", 
                       "contrast results as a ", 
                       "table."), 
                paste0("Download only the rows ", 
                       "which are left after ", 
                       "filtering (empty table if", 
                       " no filtering done)."),
                paste0("Download only the rows ", 
                       "which have been manually ", 
                       "selected (empty table if", 
                       " no selection done)."),
                paste0("This button resets any ", 
                       "rows that have been ", 
                       "selected."), 
                paste0("This checkbox shows and ", 
                       "hides options related ", 
                       "to exporting the plots."),
                paste0("Select the format for ", 
                       "exporting the plot."), 
                paste0("Specify the width of the ", 
                       "plot in pixels for export."), 
                paste0("Specify the height of the ", 
                       "plot in pixels for export."), 
                paste0("The interactive volcano ", 
                       "plot of the contrast. ", 
                       "Hovering over a data ", 
                       "point would show a ", 
                       "tooltip with additional ", 
                       "information. Clicking on ", 
                       "a data point will select ", 
                       "the respective gene in ", 
                       "the table and will also ", 
                       "generate a column plot to ", 
                       "visualise the expression ", 
                       "of the selected gene in ", 
                       "tissues. The plot can be ", 
                       "zoomed-in by clicking and ", 
                       "dragging a box over the ", 
                       "area to be zoomed. ", 
                       "Additionally the axes of ", 
                       "the plot can be adjusted ", 
                       "by dragging the edges of ", 
                       "the respective axes. The ", 
                       "zoom and/or axes can be ", 
                       "reset by double-clicking ", 
                       "in the empty space in the ", 
                       "plot area or by clicking ", 
                       "the house icon in the ", 
                       "top-right of the plot. ", 
                       "The plot can be exported ", 
                       "using the export options ", 
                       "by clicking the camera ", 
                       "icon in the top-right of ", 
                       "the plot. The same ", 
                       "functionality except for ", 
                       "data point selection also ", 
                       "exists for the column ", 
                       "plot.")
      )))
      
      observeEvent(input$TutorialButton2,{
        updateCheckboxInput(session = session, 
                            inputId = "PlotSaveOps", 
                            value = TRUE)
        introjs(session = session , options = list(steps=steps2()))
      })
    }
    
  } # closing bracket for server call
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)






