#
# A. Gaye 16th august 2016 #
# A new function to plot FDR results. This was generated using the same strategy
# as the functions 'power' and 'FDcost' in the library 'PROPER'*.
# *  Bioinformatics. 2015 Jan 15;31(2):233-41. doi: 10.1093/bioinformatics/btu640
#

plotFDR <- function (powerOutput, cols = 1:ncol(powerOutput$FDR), 
			lty = 1:ncol(powerOutput$power), 
          		main = "FDR", ylab = "FDR", leg = TRUE, error.bar = FALSE) 
{
  strata = levels(cut(0, powerOutput$strata))
  FDR.all = powerOutput$FDR
  FDR = apply(FDR.all, c(1, 2), function(x) mean(x[is.finite(x)], 
                                                 na.rm = TRUE))
  FDR.se = apply(FDR.all, c(1, 2), function(x) sd(x[is.finite(x)], 
                                                  na.rm = TRUE))/sqrt(dim(FDR.all)[3])
  plot(FDR, type = "b", lwd = 2, col = "red", lty = lty, 
          ylim = c(0, max(FDR, na.rm = TRUE)), xlim = c(0, length(FDR)+1),
          main = main, axes = FALSE, ylab = "", xlab="")
  stratas = 1:dim(powerOutput$power)[1]
  FDR.lo = FDR - FDR.se * 1.96
  FDR.hi = FDR + FDR.se * 1.96
  if (error.bar) {
    for (j in 1:1) arrows(stratas, FDR.lo, 
                                    stratas, FDR.hi, length = 0.05, angle = 90, 
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
    ypos = -diff(range(y, na.rm = TRUE)) * 0.15
    text(1:length(strata) - 0.1, ypos, strata, srt = 270 + 45, 
         xpd = TRUE, adj = c(0, 0))
  }
  add.axis1(FDR, strata)
  axis(2)
  box()
}
