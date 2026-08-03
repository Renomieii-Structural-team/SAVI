library(ggplot2)
library(dplyr)
library(colorspace)

plot_heatmap <- function(df,
                         gene,
                         save = TRUE,
                         out_dir = "objetos_ggplot",
                         palette = "Blue-Red 3",
                         dpi = 200) {
  
  size_prot <- max(df$Resno, na.rm = TRUE)
  
  # 1) define separador de eixos
  sep_scale <- case_when(
    size_prot < 600  ~ 25,
    size_prot < 1200 ~ 50,
    TRUE             ~ 100
  )
  breaks_x <- seq(1, size_prot, by = sep_scale)
  
  # 2) define dimensões do plot
  plot_dims <- switch(
    findInterval(size_prot, c(600, 1200)),
    `0` = c(6, 8),
    `1` = c(6, 12),
    `2` = c(8, 12)
  )
  
  # 3) monta o ggplot
  p <- ggplot(df, aes(x = Resno, y = ResID, fill = score)) +
    geom_tile(alpha = 0.5, color = "gray70", size = 0.1) +
    coord_fixed() +
    labs(
      title = paste0("Heatmap: ", gene),
      x     = "Amino acid index",
      y     = "Mutated Residue",
      fill  = "Score"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 8),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      panel.grid      = element_blank(),
      plot.title      = element_text(face = "bold", size = 14, hjust = 0)
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = breaks_x) +
    scale_y_discrete(expand = c(0, 0)) +
    colorspace::scale_fill_continuous_diverging(palette = palette, mid = 0.5)
  
  # 4) salva, se pedido
  if (save) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    saveRDS(p, file.path(out_dir, sprintf("heatmap_%s.rds", gene)))
    ggsave(
      filename = file.path(out_dir, sprintf("heatmap_%s.png", gene)),
      plot     = p,
      width    = plot_dims[2],
      height   = plot_dims[1],
      units    = "in",
      dpi      = dpi
    )
  }
  
  return(p)
}
