## flow repository IDs ----
flowRep_raw <- readr::read_delim("data-raw/flow_repo_links_table.txt")

flow_id <- flowRep_raw |>
  tidyr::separate(`FlowRepository ID`, sep="/", into=c(NA,NA, NA, NA, "ID")) |>
  dplyr::pull(ID)

# create hyperlinks
flow_rep_links <- paste0(
  '<a href=\"', flowRep_raw$`FlowRepository ID`, '\" target=\"_blank\">', flow_id, '</a>' 
)

flow_data_ids <- flowRep_raw |>
  dplyr::mutate(`FlowRepository ID` = flow_rep_links)

saveRDS(flow_data_ids, file="data/flow_data_ids.rds")


# Genetic Perturbation
temp_file_path <- "D:/projects/expressionViewer/www/Genetic_Perturbation/"
#temp_file_path <- "/www/Genetic_Perturbation/"
genetic_pert_names <- list.files("D:/projects/expressionViewer/www/Genetic_Perturbation/")[-8]

make_filepath_tbl <- function(genetic_pert){
  
  tissues_full <- list.files(paste0(temp_file_path, genetic_pert), full.names = TRUE)
  tissues_full <- tissues_full[!stringr::str_detect(tissues_full, "Expression_overlays")]
  all_png <- list.files(tissues_full, full.names = TRUE)
  all_png <- stringr::str_remove(all_png, "D:/projects/expressionViewer/")
  
  tibble::tibble(
    tissue = basename(tissues_full),
    tsne   = all_png[stringr::str_detect(all_png, "tsne")],
    hist   = all_png[stringr::str_detect(all_png, "histogram")]
  )
}

#make_filepath_tbl(genetic_pert_names[2])
cell_phenotype <- lapply(genetic_pert_names, make_filepath_tbl)
names(cell_phenotype) <- genetic_pert_names
genetic_perturbations <- list(cell_phenotype = cell_phenotype)


heatmap_files <- list.files(paste0(temp_file_path, "Heatmaps"), full.names = TRUE)

# CHANGE FILEPATH TO START FROM WWW


make_filepath_marker <- function(genetic_pert){
  
  markers_full <- list.files(paste0(temp_file_path, genetic_pert, "/Expression_overlays"), full.names = TRUE)
  markers <- basename(markers_full)
  markers <- stringr::str_remove(markers, pattern=".*_tsne_overlay__")
  markers <- gsub(".jpg", x = markers, replacement = "")
  

  heatmap_file <- heatmap_files[stringr::str_detect(heatmap_files, genetic_pert)]
  
  markers_full <- stringr::str_remove(markers_full, "D:/projects/expressionViewer/")
  names(markers_full) <- markers
  heatmap_file <- stringr::str_remove(heatmap_file, "D:/projects/expressionViewer/")
  
  list(markers = markers_full, heatmap = heatmap_file)
}

#make_filepath_marker(genetic_pert_names[4])
marker_files <- lapply(genetic_pert_names, make_filepath_marker)
names(marker_files) <- genetic_pert_names

genetic_perturbations[["marker_expression"]] <- marker_files

saveRDS(genetic_perturbations, file="data/genetic_perturbations.rds")
  
# Tissue expression heatmaps

#x <- list.files("D:/projects/expressionViewer_extra_data/new_data/Most highly expressed genes" , full.names = TRUE)

#sapply(x, pdftools::pdf_convert)

highly_expr_filenames <- list.files("www/Tissue_expression/highly_expr")
highly_expr_names <- gsub("_1.png", "", highly_expr_filenames)
highly_expr_names <- gsub("Top50_HEG_", "", highly_expr_names)
names(highly_expr_filenames) <- highly_expr_names
#saveRDS(highly_expr_filenames, file="data/highly_expr_filenames.rds")

# convert pdf to png
x <- list.files("D:/projects/expressionViewer_extra_data/new_data/Most differential genes/" , full.names = TRUE)

sapply(x, pdftools::pdf_convert)

diff_expr_filenames <- list.files("www/Tissue_expression/differentially_expr/")
diff_expr_names <- gsub("_vs_Blood_1.png", "", diff_expr_filenames)
diff_expr_names <- gsub("Top50_DEG_", "", diff_expr_names)
names(diff_expr_filenames) <- diff_expr_names
#saveRDS(diff_expr_filenames, file="data/diff_expr_filenames.rds")

diff_expr_filenames1 <- paste0("www/Tissue_expression/differentially_expr/", diff_expr_filenames)
names(diff_expr_filenames1) <- names(diff_expr_filenames)

highly_expr_names1 <- paste0("www/Tissue_expression/highly_expr/", highly_expr_filenames)
names(highly_expr_names1) <- names(highly_expr_filenames)


tissue_heatmap_files <- list(most_expr = as.list(highly_expr_names1), diff_expr = as.list(diff_expr_filenames1))

saveRDS(tissue_heatmap_files, file = "data/tissue_heatmap_files.rds")





