make_fig_d <- function(pal) {
  
  names(pal) <- c("Infected", "Uninfected")
  
  ann_colours <- list(Group = pal)
  
  meta <- read_delim("data/metadata.tsv") %>%
    column_to_rownames("sample_id") %>%
    mutate(Group = recode(group,
                          "Rv-Infected" = "Infected",
                          "Uninf" = "Uninfected")) %>%
    select(Group)
  
  genes <- read_delim("data/de.all_res.tsv") %>%
    rename(gene = 1) %>%
    select(gene, BH.p.value) %>%
    top_n(-50, BH.p.value) %>%
    mutate(gene = gsub("-mRNA", "", gene)) %>%
    pull(gene)
  
  mat <- read_delim("data/assay_data.tsv") %>%
    filter(gene %in% genes) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  italicised_names <- lapply(rownames(mat), function(x) bquote(italic(.(x))))
  
  title <- expression(bolditalic(z)*bold(-score))
  
  p <- ggplotify::as.ggplot(
    ComplexHeatmap::pheatmap(
      log(mat + 1),
      annotation_col = meta,
      annotation_colors = ann_colours,
      scale = "row",
      show_colnames = FALSE,
      annotation_legend = FALSE,
      labels_row = as.expression(italicised_names),
      heatmap_legend_param = list(title = title)
    )
  )
  
  return(p)
  
}
