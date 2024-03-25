get_gene_set <- function(pwy, gois) {
  
  gois %>%
    filter(pathway == pwy) %>%
    pull(gene)
  
}

run_gsea <- function(de, gene_set) {
  
  dat <- read_delim(de,
                    col_select = 1:2) %>%
    setNames(c("hgnc_symbol", "logFC")) %>%
    mutate(hgnc_symbol = gsub("-mRNA$", "", hgnc_symbol))
  
  ranks <- dat$logFC
  names(ranks) <- dat$hgnc_symbol
  
  set.seed(123456)
  res <- fgsea::fgsea(gene_set, ranks)
  
  return(res)
  
}

make_fig_e <- function(pal) {
  
  gois <- read_delim("data/gois.txt") %>%
    filter(!grepl("Complex", pathway) & pathway != "Fatty Acid Metabolism")
  
  pathways <- unique(gois$pathway)
  
  gene_set <- pathways %>%
    purrr::map(\(x) get_gene_set(x, gois)) %>%
    setNames(pathways)
  
  res <- run_gsea("data/de.all_res.tsv", gene_set) %>%
    mutate(group = ifelse(NES >= 0,
                          paste("Infected"),
                          paste("Uninfected"))) %>%
    mutate(label = ifelse(padj <= 0.05, "*", "")) %>%
    mutate(hjust = ifelse(NES > 0, -0.25, 1.25))
  
  pal_names <- c("Infected", "Uninfected")
  
  names(pal) <- pal_names
  
  p <- ggplot(res, aes(x = NES, y = reorder(pathway, NES))) +
    geom_bar(aes(fill = group),
             stat = "identity",
             colour = "black",
             width = 1) +
    geom_text(aes(label = label,
                  hjust = hjust,
                  size = -log10(pval)),
              show.legend = FALSE,
              vjust = 0.75) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = pal) +
    labs(x = "Normalised enrichment score", y = "Pathway", fill = "") +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.text = ggtext::element_markdown(),
          legend.position = "left",
          legend.title = element_blank()) +
    scale_size_continuous(range = c(3, 10)) +
    scale_x_continuous(expand = expansion(mult = 0.15, 0.15))
  
  return(p)
  
}
