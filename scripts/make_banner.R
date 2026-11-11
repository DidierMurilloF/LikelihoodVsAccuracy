# Banner for the article: a stylized REML "likelihood hill".
# Draws the profile REML log-likelihood of the example dataset
# (results/curves_iid.rds) against heritability and writes a 1600 x 600 px
# image to images/banner.png, plus a 1200 x 630 px social preview (the
# `preview:` image of index.Rmd) to images/preview.png.
# From the repo root:   Rscript scripts/make_banner.R

suppressPackageStartupMessages(library(ggplot2))

curve <- readRDS("results/curves_iid.rds")
gamma_hat <- attr(curve, "gamma_hat")
h2_hat <- gamma_hat / (1 + gamma_hat)

# Keep the top of the hill (within 12 log-likelihood units of the peak) and
# rescale its height to [0, 1].
hill <- subset(curve, loglik >= max(loglik) - 12)
hill$height <- (hill$loglik - min(hill$loglik)) / diff(range(hill$loglik))
peak <- approx(hill$h2, hill$height, xout = h2_hat)$y

bg <- "#0f1b2d"
ink <- "#f0efec"
accent <- "#eb6834"
hill_col <- "#3987e5"

# The same picture at any size; title and subtitle sizes are in mm (ggplot2
# text units), so they are set per image.
hill_plot <- function(title_size, subtitle_size, title_y, subtitle_y, margin) {
  ggplot(hill, aes(h2, height)) +
    geom_area(fill = hill_col, alpha = 0.25) +
    geom_line(colour = hill_col, linewidth = 1.6) +
    annotate("segment", x = h2_hat, xend = h2_hat, y = title_y + 0.12, yend = peak,
             colour = accent, linewidth = 0.9, linetype = "22") +
    annotate("point", x = h2_hat, y = peak, colour = accent, size = 5) +
    annotate("text", x = 0.02, y = title_y, colour = ink,
             label = "Likelihood vs Prediction Accuracy",
             size = title_size, hjust = 0, fontface = "bold") +
    annotate("text", x = 0.02, y = subtitle_y, hjust = 0, colour = "#c3c2b7",
             label = "A higher peak is not a better prediction",
             size = subtitle_size) +
    scale_x_continuous(limits = c(0, max(hill$h2)), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = bg, colour = NA),
          panel.background = element_rect(fill = bg, colour = NA),
          plot.margin = margin)
}

dir.create("images", showWarnings = FALSE)
banner <- hill_plot(17, 9, 0.30, 0.14, margin(30, 40, 0, 40))
ggsave("images/banner.png", banner, width = 16, height = 6, dpi = 100, units = "in")
cat("Wrote images/banner.png\n")

preview <- hill_plot(13.5, 7.5, 0.30, 0.15, margin(30, 40, 0, 40))
ggsave("images/preview.png", preview, width = 12, height = 6.3, dpi = 100, units = "in")
cat("Wrote images/preview.png\n")
