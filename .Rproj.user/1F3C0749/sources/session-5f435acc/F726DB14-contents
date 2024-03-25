imp_dat <- function(file) {
  
  read_tsv(file) %>%
    pivot_longer(!gene, names_to = "sample_id", values_to = "expr") %>%
    select(1:3) %>%
    mutate(sample_id = str_pad(sample_id, 2, pad = "0"))
  
}

import_data <- function() {
  
  metadata <- read_delim("data/metadata.tsv") %>%
    mutate(sample_id = str_pad(sample_id, 2, pad = "0")) %>%
    rename(disease = group) %>%
    mutate(disease = recode(disease,
                            "Rv-Infected" = "Infected",
                            "Uninf" = "Uninfected"))
  
  dat <- read_delim("data/assay_data.tsv") %>%
    pivot_longer(!gene, names_to = "sample_id", values_to = "expr") %>%
    mutate(sample_id = str_pad(sample_id, 2, pad = "0")) %>%
    inner_join(metadata, by = "sample_id")
  
  return(dat)
  
}

make_boxplot <- function(ldh, GROUP, method, method.args = list(), pal) {
  
  if (GROUP == "disease") {
    
    p_inp <- ldh %>%
      select(disease, gene, expr) %>%
      setNames(c("x", "y", "expr")) %>%
      mutate(x = recode(x, "Infected" = "Inf.", "Uninfected" = "Uninf.")) %>%
      mutate(x = factor(x, levels = c("Uninf.", "Inf.")))
    
    strip_text <- "bold.italic"
    
    axis_text_x <- element_text()
    
    title <- NULL
    
    names(pal) <- c("Inf.", "Uninf.")
    
  } else {
    
    p_inp <- ldh %>%
      select(gene, disease, expr) %>%
      setNames(c("x", "y", "expr"))
    
    strip_text <- "bold"
    
    axis_text_x <- element_text(face = "italic")
    
  }
  
  ggplot(p_inp, aes(x = x, y = log(expr))) +
    facet_grid(~ y, scales = "free") +
    geom_boxplot(aes(colour = x),
                 show.legend = FALSE,
                 outlier.shape = NA) +
    geom_jitter(aes(colour = x, fill = x),
                pch = 21,
                alpha = 0.5,
                show.legend = FALSE) + 
    theme_bw(base_size = 12.5) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          axis.text.x = axis_text_x,
          strip.text = element_text(face = strip_text),
          axis.title = element_text(face = "bold")) +
    labs(x = "Gene", y = "log(count)") +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    ggpubr::stat_pwc(method = method,
                     method.args = method.args,
                     p.adjust.method = "BH") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.105)))
  
}

plot_dat <- function(method, method.args = list(), pal) {
  
  dat <- import_data()
  
  ldh <- dat %>%
    filter(gene %in% c("LDHA", "LDHB")) %>%
    select(disease, gene, expr) %>%
    distinct()
  
  ldh$disease <- factor(ldh$disease, levels = c("Uninfected", "Infected"))
  
  p1 <- make_boxplot(ldh, "disease", method, method.args, pal = pal)
  
  p2 <- make_boxplot(ldh,
                     "gene",
                     method,
                     method.args,
                     pal = wesanderson::wes_palette("Chevalier1", 2))
  
  p <- patchwork::wrap_plots(p1, p2)
  
  return(p)
  
}

make_fig_g <- function(pal) {
  plot_dat(method = "wilcox_test", method.args = list(paired = TRUE), pal) 
}
