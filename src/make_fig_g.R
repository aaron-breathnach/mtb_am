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
    labs(x = "Condition", y = "log(count)") +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    ggpubr::stat_pwc(method = method,
                     method.args = method.args,
                     p.adjust.method = "BH") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.105)))
  
}

make_bar_chart <- function(ldh, GROUP, method, method.args = list()) {
  
  df_1 <- ldh %>%
    mutate(expr = log1p(expr)) %>%
    mutate(disease = factor(disease, levels = c("Uninfected", "Infected")))
  
  df_2 <- df_1 %>%
    group_by(disease, gene) %>%
    summarise(mean = mean(expr), se = sd(expr) / sqrt(n()))
  
  ann <- df_1 %>%
    group_by(disease) %>%
    rstatix::wilcox_test(expr ~ gene, paired = TRUE) %>%
    mutate(gene = group1)
  
  y <- max(df_1$expr) + (0.1 * ((1.1 * max(df_1$expr)) - (min(df_1$expr) / 1.1)))
  
  ggplot(df_2, aes(x = gene, y = mean)) +
    facet_wrap(~ disease) +
    geom_bar(
      stat = "identity",
      colour = "black",
      fill = "#967bb6"
    ) +
    geom_jitter(
      data = df_1,
      aes(x = gene, y = expr),
      colour = "grey",
      alpha = 0.500,
      width = 0.125
    ) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      width = 0.25
    ) +
    theme_bw(base_size = 12.5) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "italic"),
      strip.text = element_text(face = "bold")
    ) +
    coord_cartesian(ylim = c(min(df_1$expr), 1.05 * max(df_1$expr))) +
    labs(x = "Gene", y = "log(count)") +
    ggpubr::stat_pvalue_manual(
      ann,
      label = "p = {round(p, 3)}",
      y.position = y
    )
  
}

plot_dat <- function(method, method.args = list(), pal) {
  
  dat <- import_data()
  
  ldh <- dat %>%
    filter(gene %in% c("LDHA", "LDHB")) %>%
    select(disease, gene, expr) %>%
    distinct()
  
  ldh$disease <- factor(ldh$disease, levels = c("Uninfected", "Infected"))
  
  p1 <- make_boxplot(ldh, "disease", method, method.args, pal = pal)
  
  p2 <- make_bar_chart(ldh, "disease", method, method.args)
  
  p <- p1 + p2
  
  return(p)
  
}

plot_dat <- function(method, method.args = list(), pal) {
  
  dat <- import_data()
  
  ldh <- dat %>%
    filter(gene %in% c("LDHA", "LDHB")) %>%
    select(disease, gene, expr) %>%
    distinct()
  
  ldh$disease <- factor(ldh$disease, levels = c("Uninfected", "Infected"))
  
  p1 <- make_boxplot(ldh, "disease", method, method.args, pal = pal)
  
  p2 <- make_bar_chart(ldh, "disease", method, method.args)
  
  p <- patchwork::wrap_plots(p1, p2)
  
  return(p)
  
}

make_fig_g <- function(pal) {
  plot_dat(method = "wilcox_test", method.args = list(paired = TRUE), pal) 
}
