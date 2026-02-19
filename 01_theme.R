library(ggplot2)
library(scales)

color_gam     <- "#16A085"
color_classic <- "#C0392B"

method_colors <- c(
  "GAM" = color_gam, "GAM-logIOR" = color_gam, "GAM-RERI" = color_gam, 
  "GAM-Doble" = color_gam, "GAM-IOR" = color_gam,
  "Estratificado" = color_classic, "Estratificado-IOR" = color_classic, 
  "Estratificado-RERI" = color_classic, "Estratificado-Doble" = color_classic
)

theme_base <- function() {
  theme_bw(base_size = 15) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "gray90", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold", size = 15, margin = margin(3, 5, 3, 5)),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.8, "lines"),
      text = element_text(face = "bold"),
      axis.text = element_text(face = "bold", size = 15),
      axis.title = element_text(face = "bold", size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )
}

theme_set(theme_base())

scale_method <- function(aesthetics = c("color", "fill")) {
  aesthetics <- match.arg(aesthetics)
  if (aesthetics == "color") scale_color_manual(values = method_colors)
  else scale_fill_manual(values = method_colors)
}