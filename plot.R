# Define custom function to plot the section
plot_section <- function() {
  # Plotting the section between [0, 2000]
  x1 <- seq(0, 2000, length.out = 100)
  y1 <- rep(1, length(x1))
  polygon(x = c(0, x1, 2000), y = c(0, y1, 0), col = "blue", border = NA)
  
  # Plotting the section between [2000, 10000]
  x2 <- seq(2000, 10000, length.out = 100)
  y2 <- rep(1, length(x2))
  polygon(x = c(2000, x2, 10000), y = c(0, y2, 0), col = "red", border = NA)
  
  # Add title and axis labels
  title(main = "Section Plot", xlab = "X", ylab = "Y")
}

# Plot the section
plot_section()
