gene_list <- function(pwy, gois) {
  gois %>%
    filter(pathway == pwy) %>%
    pull(gene)
}

viz_log_fc <- function(pwy, gois, de) {
  
  genes <- gene_list(pwy, gois)
  
  df <- de %>%
    filter(gene %in% genes)
  
  n <- 0.05 * abs(max(df$log_fc) - min(df$log_fc))
  
  df <- df %>%
    mutate(sig_pos = ifelse(log_fc < 0, log_fc - n, log_fc + n))
  
  val <- ifelse(min(df$log_fc) < 0, 0.1, 0)
  
  n_feat <- nrow(df)
  row_height <- n_feat * 0.2
  fig_height <- n_feat * 0.2 + 2
  
  p <- ggplot(df, aes(x = log_fc, y = reorder(gene, log_fc))) +
    geom_bar(aes(fill = -log10(p_value)), stat = "identity", width = 1, colour = "black") +
    geom_text(aes(x = sig_pos, y = reorder(gene, log_fc), label = sig),
              vjust = 0.75,
              size = 7.5) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme_classic(base_size = 12.5) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title.x = ggtext::element_markdown(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          axis.text.y = element_text(face = "italic"),
          legend.title = ggtext::element_markdown(face = "bold"),
          legend.position = "top") +
    labs(title = pwy,
         x = "log<sub>2</sub> fold change",
         y = "Gene",
         fill = "-log<sub>10</sub>(*p*-value)") +
    scale_x_continuous(expand = expansion(mult = c(val, 0.1))) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    ggh4x::force_panelsizes(rows = unit(nrow(df) * 0.2, "in"),
                            cols = unit(3, "in"))
  
  return(p)
  
}

make_fig_f <- function() {
  
  pathways <- c("Glycolysis", "TCA Cycle", "Mitochondrial OXPHOS")
  
  gois <- read_delim("data/gois.txt") %>%
    filter(pathway %in% pathways)
  
  de <- read_delim("data/de.all_res.tsv") %>%
    select(1, 2, 9, 10, 12) %>%
    setNames(c("gene", "log_fc", "p_value", "fdr", "gene_set")) %>%
    mutate(gene = gsub("-mRNA", "", gene)) %>%
    mutate(sig = ifelse(fdr <= 0.05, "*", ""))
  
  plot_list <- lapply(pathways, viz_log_fc, gois, de)
  
  p <- patchwork::wrap_plots(plot_list)
  
  return(p)
  
}

