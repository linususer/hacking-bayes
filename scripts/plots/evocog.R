plot_spaced_barplot <- function(names, values, title, spacing = 2) {
  if (length(names) != length(values)) {
    stop("Length of 'names' and 'values' must be equal.")
  }
  
  # Prepare spaced vectors
  values_spaced <- c()
  names_spaced <- c()
  
  for (i in seq_along(values)) {
    values_spaced <- c(values_spaced, values[i])
    names_spaced <- c(names_spaced, names[i])
    
    if (i %% spacing == 0 && i != length(values)) {
      values_spaced <- c(values_spaced, 0)      # dummy bar
      names_spaced <- c(names_spaced, "")       # empty label
    }
  }
  
  # Plot with dummy bars hidden (NA color)
  barplot(values_spaced,
          names.arg = names_spaced,
          col = ifelse(values_spaced == 0, NA, "steelblue"),
          border = NA,
          main = title,
          xlab = "",
          ylab = "Reaction Time (RT)",
          axes = FALSE,
          width = 0.25,
          space = 1.2
          )
}

# First Block
fnames <- c("Baseline", "TBR", "TBF")
fvalues <- c(3, 1, 4)
pdf("evocog_barplot.pdf")
plot_spaced_barplot(fnames, fvalues, "First Part: Experiment by Dames & Oberauer (2022)", spacing = 3)
dev.off()
# Second Block
fnames2 <- c("Baseline", "TBR (vis[vis])", "TBR (aud[aud])", "TBR (vis[aud])", "TBR (aud[vis])", "TBF (vis)", "TBF (aud)")
fvalues2 <- c(3,1,1,3,3,4,4)
pdf("evocog_barplot2.pdf", width = 10, height = 6)
plot_spaced_barplot(fnames2, fvalues2, "Second Part: Multimodal Extension", spacing = 7)
dev.off()
fvalues2 <- c(3,1,1,1,1,4,4)
pdf("evocog_barplot3.pdf", width = 10, height = 6)
plot_spaced_barplot(fnames2, fvalues2, "Second Part: Multimodal Extension", spacing = 7)
dev.off()
