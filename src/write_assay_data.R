get_sample_id <- function(rcc) {
  
  tmp <- rcc %>%
    str_split("_") %>%
    unlist() %>%
    nth(3) 
  
  str_c(
    str_sub(tmp, 1, 5),
    "_",
    ifelse(grepl("Un", tmp), "uninfected", "infected")
  )
  
}

import_assay_data <- function(inp_dir) {
  
  inpdir <- gsub("/$", "", inp_dir)
  
  rccs <- list.files(inpdir,
                     pattern = ".RCC",
                     full.names = TRUE,
                     recursive = TRUE)
  
  rcc_obj <- NanoStringNCTools::readNanoStringRccSet(rccs)
  
  genes <- Biobase::fData(rcc_obj) %>%
    pull(GeneName)
  
  assay_data <- Biobase::assayDataElement(rcc_obj, "exprs")
  
  rownames(assay_data) <- genes
  
  sample_ids <- colnames(assay_data) %>%
    purrr::map(function(x) get_sample_id(x)) %>%
    unlist()
  
  colnames(assay_data) <- sample_ids
  
  assay_data <- assay_data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    as_tibble()
  
  return(assay_data)
  
}

write_assay_data <- function() {
  
  assay_data <- import_assay_data("data/rcc")
  
  write_tsv(assay_data, "data/assay_data.tsv")
  
}
