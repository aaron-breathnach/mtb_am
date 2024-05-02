library(tidyverse)

source("src/write_assay_data.R")

############
## SERIES ##
############

series <- read_tsv("data/geo_sub/series.txt")

#############
## SAMPLES ##
#############

files <- list.files("data/rcc", full.names = TRUE)

raw_data_file <- basename(files)

sample_name <- purrr::map(files, function(x) get_sample_id(x)) %>%
  unlist()

tab <- tibble("Sample name" = sample_name, "raw data file" = raw_data_file)

meta <- read_delim("data/metadata.tsv") %>%
  mutate(group = ifelse(group == "Uninf", "uninfected", "infected")) %>%
  mutate(title = sprintf("%s %s", subject_id, group)) %>%
  mutate(smoking = ifelse(smoking == "N", "never-smoker", "ex-smoker")) %>%
  dplyr::rename(`Sample name` = 1) %>%
  mutate(`source name` = "Human alveolar macrophages",
         organism = "Homo sapiens",
         molecule = "RNA",
         label = "n/a",
         description = series[[3, 2]],
         platform = "nCounter Metabolic Pathways Panel") %>%
  inner_join(tab, by = "Sample name")

cols <- c("subject_id", "group", "smoking")

characteristics <- cols %>%
  purrr::map(function(x) paste("characteristics:", str_replace(x, "_", " "))) %>%
  unlist()

meta_1 <- meta %>%
  select(`Sample name`, title, `raw data file`, `source name`, organism)

meta_2 <- meta %>%
  select(subject_id, group, smoking) %>%
  setNames(characteristics)

meta_3 <- meta %>%
  select(molecule, label, description, platform)
  
samples <- purrr::reduce(list(meta_1, meta_2, meta_3), cbind) %>%
  as_tibble()

###############
## PROTOCOLS ##
###############

protocols <- read_tsv("data/geo_sub/protocols.txt", col_names = FALSE)

##################
## save to xlsx ##
##################

wb <- openxlsx::createWorkbook()

## sheet 1

openxlsx::addWorksheet(wb, sheetName = "Metadata")

title_1 <- tibble(X1 = "SERIES")
title_2 <- tibble(X2 = "SAMPLES")
title_3 <- tibble(X3 = "PROTOCOLS")

l <- list(title_1, series, title_2, samples, title_3, protocols)

start_row <- c(
  1, # title_1
  2, # series
  2 + nrow(series) + 1, # title_2
  2 + nrow(series) + 1 + 1, # samples
  2 + nrow(series) + 1 + 1 + nrow(samples) + 2, # title_3
  2 + nrow(series) + 1 + 1 + nrow(samples) + 2 + 1 # protocols
)

col_names <- c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)

for (i in 1:length(l)) {
  
  openxlsx::writeData(
    wb,
    sheet = "Metadata",
    x = l[[i]],
    startCol = 1,
    startRow = start_row[i],
    colNames = col_names[i]
  )
  
}

## sheet 2

Matrix <- read_delim("data/assay_data.tsv")

openxlsx::addWorksheet(wb, sheetName = "Matrix")

openxlsx::writeData(
  wb,
  sheet = "Matrix",
  x = Matrix
)

out_dir <- "data/geo_sub/submitted"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

xlsx <- sprintf("%s/geo_sub.xlsx", out_dir)
openxlsx::saveWorkbook(wb, xlsx, overwrite = TRUE)

#####################
## write raw files ##
#####################

cp <- function(rcc) {
  cmd <- sprintf("cp %s data/geo_sub/submitted/", rcc)
  system(cmd)
}

rcc <- Sys.glob("data/rcc/*.RCC")

lapply(rcc, cp)

system("cd data/geo_sub; zip submitted.zip submitted/*; rm -r submitted")

