make_figure_1 <- function(format) {
  
  fig_a <- make_fig_a()
  fig_b <- make_fig_b(pal = c("Infected" = "#E7B5AC", "Uninfected" = "#869F77"))
  fig_c <- make_fig_c(pal = c("red", "blue"))
  fig_d <- make_fig_d(pal = c("#E7B5AC", "#869F77"))
  fig_e <- make_fig_e(pal = c("red", "blue"))
  fig_f <- make_fig_f()
  fig_g <- make_fig_g(pal = c("#E7B5AC", "#869F77"))
  
  dir.create("figure_1", showWarnings = FALSE)
  ggsave("figure_1/fig_b.png", fig_b, width = 6.25, height = 5.00)
  ggsave("figure_1/fig_c.png", fig_c, width = 6.25, height = 5.00)
  ggsave("figure_1/fig_d.png", fig_d, width = 6.25, height = 10.0)
  ggsave("figure_1/fig_e.png", fig_e, width = 8.75, height = 5.00)
  ggsave("figure_1/fig_f.png", fig_f, width = 12.5, height = 7.50)
  ggsave("figure_1/fig_g.png", fig_g, width = 8.75, height = 5.00)
  
  if (format == "wide") {
    
    part_1 <- cowplot::plot_grid(plot.new(), fig_b, fig_c, plot.new(),
                                 nrow = 4,
                                 rel_heights = c(0.2, 1, 1, 0.2),
                                 align = "v",
                                 axis = "lr",
                                 scale = 0.95)
    
  } else {
    
    part_1 <- cowplot::plot_grid(fig_b, fig_c,
                                 nrow = 2,
                                 align = "v",
                                 axis = "lr",
                                 scale = 0.95)
    
  }
  
  part_2 <- fig_d
  
  part_1_2 <- cowplot::plot_grid(fig_a, part_1, part_2, 
                                 nrow = 1,
                                 rel_widths = c(1.5, 1, 1),
                                 scale = 0.95)
  
  e_g <- cowplot::plot_grid(fig_e, fig_g,
                            nrow = 1,
                            align = "h",
                            axis = "tb",
                            rel_widths = c(0.875, 1),
                            scale = 0.95)
  
  part_3 <- cowplot::plot_grid(e_g, fig_f,
                               nrow = 2,
                               rel_heights = c(0.5, 1),
                               scale = 0.95)
  
  if (format == "wide") {
    
    fig_1 <- cowplot::plot_grid(part_1_2, part_3,
                                nrow = 1,
                                rel_widths = c(1, 1),
                                scale = 0.95)
    
    w <- 35.0
    h <- 12.5
    
  } else {
    
    fig_1 <- cowplot::plot_grid(part_1_2, part_3,
                                nrow = 2,
                                rel_heights = c(1, 1.5),
                                scale = 0.95)
    
    w <- 17.5
    h <- 20.0
    
  }
  
  filename <- sprintf("figure_1/fig_1_%s.png", format)
  
  ggsave(filename, fig_1, width = w, height = h, bg = "white")
  
}
