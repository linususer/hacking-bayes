library(BayesFactor)

point_and_point_prior <- function() {
  # plot point prior^2
  pdf("figures/report/point-prior-and-point-prior.pdf", width = 8, height = 4)
  # double the font size
  plot(0, 0,
    xlim = c(0, 1), ylim = c(0, 1), type = "n",
    main = "Point Prior vs Point Prior",
    ylab = "", xlab = bquote("Probability of the coin for heads " * p),
    cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5
  )

  # Plot the point prior densities
  points(0.5, 0.5, col = "#CD1076", pch = 19, cex = 2)
  lines(c(0.5, 0.5), y = c(0, 0.5), col = "#CD1076", lwd = 4)
  points(0.75, 0.5, col = "#59B3E6", pch = 19, cex = 2)
  lines(c(0.75, 0.75), y = c(0, 0.5), col = "#59B3E6", lwd = 4)
  legend("topleft", legend = c(bquote(H[0]:p == 0.5), bquote(H[1]:p == 0.75)), fill = c("#CD1076", "#59B3E6"))
  dev.off()
}
point_and_point_prior()

fixed_opt_random_walk <- function() {
  set.seed(0)
  n_max <- 50
  effect_size <- 0.2
  x <- c(rnorm(1, effect_size, 1))
  bf_crit <- 3
  r_scale <- 1 / sqrt(2)
  count <- 2
  opt_done <- FALSE
  bf_opt <- c()
  bf_fixed <- c()
  res_opt <- c()

  while (count <= n_max) {
    x <- c(x, rnorm(1, effect_size, 1))
    bf_fixed <- c(bf_fixed, 1 / extractBF(ttestBF(x, mu = 0, r = r_scale))$bf)
    bf_opt <- c(bf_opt, 1 / extractBF(ttestBF(x, mu = 0, r = r_scale))$bf)
    count <- count + 1
    if (opt_done) {
      next
    } else if (tail(bf_opt, n=1L) > bf_crit) {
      res_opt <- bf_opt
      opt_done <- TRUE
    } else if (tail(bf_opt, n = 1L) < 1 / bf_crit) {
      res_opt <- bf_opt
      opt_done <- TRUE
    }
  }
  print(paste(length(x), length(bf_fixed), length(bf_opt), length(res_opt), sep=" "))
  print(tail(bf_fixed))
  print(tail(bf_opt))
  pdf("figures/report/fixed-and-optional-stopping-random-walk.pdf", height = 8, width = 12)
  par(mar = c(5, 5, 5, 5))
  plot(2:n_max, bf_fixed, type = "l", xlab = "number of trials n", ylab = bquote("Bayes Factor " * BF["01"]), col = "orange",
  ylim = c(0,3.5), xlim = c(0,50), lwd = 5, lty = 4, cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5,
  main = "Fixed Sample Size vs Optional Stopping")
  lines(2:(length(res_opt)+1), res_opt, col = "black", lwd = 5, lty = 1)
  abline(h=bf_crit, lty=2, lwd = 5)
  abline(h=1/bf_crit, lty = 2, lwd = 5)
  points(length(res_opt) + 1, tail(res_opt, n=1L), col = "black", pch=4, lwd = 5, cex = 3)
  points(length(bf_fixed) + 1, tail(bf_fixed, n=1L), col = "orange", pch=4, lwd = 5, cex = 3)
  legend("topright", legend = c("optional stopping", "fixed sample size"), lty=c(1,4), col=c("black", "orange"), lwd=5, cex=1.5)
  dev.off()
}
fixed_opt_random_walk()
