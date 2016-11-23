#
# A. Gaye 16th June 2015 #
# An amended version of the function 'plotPower' to allow for the plotting
# of just one curve. The original function throws an 'xlim error' when using 
# one sample size rather than two or more as set by he author. The plot was
# also amended to enable line and dots along with other minor changes.
#

plotPower_ag <- function (powerOutput, cols = 1:ncol(powerOutput$FD), lty = 1:ncol(powerOutput$power), 
          main = "Power", ylab = "Power", leg = TRUE, error.bar = FALSE) 
{
  nsims = dim(powerOutput$power)[3]
  power = apply(powerOutput$power, c(1, 2), mean, na.rm = TRUE)
  power.se = apply(powerOutput$power, c(1, 2), sd, na.rm = TRUE)/sqrt(nsims)
  ix.na = apply(power, 1, function(x) all(is.na(x)))
  power = power[!ix.na, ]
  power.se = power.se[!ix.na, ]
  strata = levels(cut(0, powerOutput$strata))
  strata = strata[!ix.na]

  plot(power, type = "b", lwd = 2, col = "red", lty = lty, 
          ylim = c(0, max(power, na.rm = TRUE)), xlim = c(0, length(power)+1),
          main = main, axes = FALSE, ylab = "", xlab="")
  stratas = 1:dim(powerOutput$power)[1]
  power.lo = power - power.se * 1.96
  power.hi = power + power.se * 1.96
  if (error.bar) {
    for (j in 1:1) arrows(stratas, power.lo, 
                                    stratas, power.hi, length = 0.05, angle = 90, 
                                    code = 3, col = j, lty = 1, lwd = 1)
  }
  xlab <- switch(powerOutput$stratify.by, expr = "Average count strata", 
                 dispersion = "Dispersion strata")
  mtext(xlab, side = 1, line = 4)
  mtext(ylab, side = 2, line = 2.5)
  if (leg) 
    legend("bottomright", legend = paste("n =", powerOutput$Nreps), 
           col="red", lty = c(lty), lwd = 2, bty='n')
  grid()
  add.axis1 <- function (y, strata) 
  {
    tmp = axis(1, 1:length(strata), labels = FALSE)
    ypos = -diff(range(y, na.rm = TRUE)) * 0.10
    text(1:length(strata) - 0.1, ypos, strata, srt = 270 + 45, 
         xpd = TRUE, adj = c(0, 0))
  }
  add.axis1(power, strata)
  axis(2)
  box()

}
