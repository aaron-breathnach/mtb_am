pca_plot <- function(pca_inp, meta, col_by, pal, conf = FALSE) {
  
  pca <- prcomp(pca_inp, scale. = TRUE)
  
  pca_summary <- summary(pca)$importance[2,]
  pc1_label <- paste0("PC1 (", 100 * round(pca_summary[[1]], 2), "%)")
  pc2_label <- paste0("PC2 (", 100 * round(pca_summary[[2]], 2), "%)")
  
  pca_points <- pca$x[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    as_tibble() %>%
    mutate(sample_id = as.character(sample_id)) %>%
    inner_join(meta, by = "sample_id")
  
  pca_plot <- ggplot(pca_points, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = !!ensym(col_by)), pch = 21, colour = "black", size = 3) +
    geom_point(pch = 21, size = 3) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.position = "top") +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    labs(x = pc1_label, y = pc2_label, colour = "Condition", fill = "Condition") +
    stat_ellipse(aes(colour = !!ensym(col_by)), linetype = "dashed")
  
  return(pca_plot)
  
}

pca_wrapper <- function(inp_dir, pal, conf = FALSE) {
  
  inp_dir <- paste0(gsub("/$", "", inp_dir), "/")
  
  meta <- read_tsv(paste0(inp_dir, "metadata.tsv")) %>%
    mutate(sample_id = str_pad(sample_id, 2, "left", "0")) %>%
    mutate(group = recode(group, "Uninf" = "Uninfected", "Rv-Infected" = "Infected")) %>%
    mutate(group = factor(group, levels = c("Uninfected", "Infected")))
  
  pca_inp <- read_tsv(paste0(inp_dir, "assay_data.tsv")) %>%
    pivot_longer(!gene, names_to = "sample_id", values_to = "count") %>%
    pivot_wider(names_from = "gene", values_from = "count") %>%
    column_to_rownames("sample_id")
  
  p <- pca_plot(pca_inp, meta, col_by = "group", pal)
  
  return(p)
  
}

make_fig_b <- function(pal) {
  
  pca_wrapper("data", pal, conf = TRUE)
  
}
