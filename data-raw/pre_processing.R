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





