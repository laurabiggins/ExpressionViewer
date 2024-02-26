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

# FlowRepID <- c(
#   '<a href="http://flowrepository.org/id/FR-FCM-Z6L9" target="_blank">FR-FCM-Z6L9</a>',
#   '<a href="http://flowrepository.org/id/FR-FCM-Z6LT" target="_blank">FR-FCM-Z6LT</a>'
# )

# flow_data_ids <- tibble::tibble(
#   Dataset = c("Ageing", "Microbiome"), 
#   `FlowRepository ID` = FlowRepID
# )

# DT::datatable(
#   flow_data_ids,
#   escape   = FALSE,
#   rownames = FALSE, 
#   options  = list(dom = "t")
# )

