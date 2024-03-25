import_DE <- function(DE) {
  
  cols <- c("Log2 fold change", "P-value", "BH.p.value")
  
  tab <- read_delim(DE) %>%
    select(1, all_of(cols)) %>%
    setNames(c("gene", "log2FoldChange", "pvalue", "padj")) %>%
    mutate(gene = gsub("-.*", "", gene))
  
  return(tab)
  
}

volcano_plot <- function(de, q = 0.05, pal) {
  
  groups <- c(
    paste("\U2191", "in infected AM"),
    paste("\U2193", "in infected AM"),
    "non-significant"
  )
  
  tab <- de %>%
    mutate(colour = case_when(
      padj <  0.05 & log2FoldChange > +0001 ~ groups[1],
      padj <  0.05 & log2FoldChange < -0001 ~ groups[2],
      padj > 0.05 | abs(log2FoldChange) < 1 ~ groups[3]
    )) %>%
    mutate(colour = factor(colour, levels = groups))
  
  x_min <- min(tab$log2FoldChange)
  x_max <- max(tab$log2FoldChange)
  x_val <- max(abs(c(x_min, x_max)))
  x_lim <- c(-x_val, x_val)
  
  ann <- tab %>%
    drop_na() %>%
    filter(abs(log2FoldChange) > 2) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "up", "down")) %>%
    top_n(-10, padj)
  
  names(pal) <- groups
  
  tmp_max_sig <- tab %>%
    filter(padj <= q)
  
  tmp_volcano <- ggplot(tab, aes(x = log2FoldChange, y = -log10(pvalue)))
  
  if (nrow(tmp_max_sig) > q) {
    
    max_sig <- filter(tmp_max_sig, pvalue == max(pvalue)) %>%
      .[[1, "pvalue"]]
    
    tmp_volcano <- tmp_volcano +
      geom_hline(yintercept = -log10(max_sig), linetype = "dashed", colour = "darkgrey")
    
  }
  
  volcano <- tmp_volcano +
    geom_vline(xintercept = c(-2, -1, 1, 2), linetype = "dashed", colour = "darkgrey") +
    geom_point(aes(colour = colour)) +
    theme_bw(base_size = 12.5) +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = ggtext::element_markdown(),
          axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          panel.border = element_blank(),
          legend.position = "top") +
    scale_colour_manual(values = pal, breaks = groups[1:2]) +
    labs(x = "**log<sub>2</sub> fold change**",
         y = "**-log<sub>10</sub>(*p*)**",
         colour = NULL) +
    xlim(x_lim) +
    ggrepel::geom_text_repel(data = ann,
                             aes(x = log2FoldChange, y = -log10(pvalue), label = gene),
                             fontface = "italic",
                             min.segment.length = 0) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  
  return(volcano)
  
}

make_fig_c <- function(pal) {
  
  pal <- purrr::map(pal, function(x) colorspace::darken(x, 1/3)) %>%
    unlist()
  
  pal <- c(pal, "lightgrey")
  
  de <- import_DE("data/de.all_res.tsv")
  
  p <- volcano_plot(de, q = 0.05, pal)
  
  return(p)
  
}
