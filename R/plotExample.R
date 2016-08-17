#
# A. Gaye 16th August 2016 #
# Example on how to use the new and amended functions 
# to compute and plot the power from the information
# obtained from the simulation (see script 'simExample.R').
# N.B: This is just an example and the results are just an 
# illustration have no real significance.
#
  
# load the library
library("PROPER")

# load the object that holds the RNA-seq DE detection simulation rersults
load("runSimsResults.Rdata")

# estimate marginal power without exclusions
powers <- comparePower(runSimsResults, alpha.type="fdr", alpha.nominal=0.05, 
                       stratify.by="expr", delta=0.585, filter.by="none", 
                       strata=c(0, 5, 10, 2^(1:6) * 10, Inf), strata.filtered=1)
# plot power by strata: BE AWARE THAT WHEN STRATA IS EXCLUDED 
# THE CURVE SHOULD SHIFT RIGHT BUT IS DOES NOT - THIS NEEDS ADDRESSING

# load the two functions to use
source("plotPower_ag.R")
source("plotFDR.R")
par(mfcol=c(1,2), oma = c(0, 0, 2, 0))
toptitle <- "OUTCOME"
plotPower_ag(powers, main=toptitle, error.bar=TRUE)
plotFDR(powers, main=toptitle, error.bar=TRUE)

# store details
p0 <- summaryPower(powers)

# estimate marginal power after excluding genes with mean read counts <= 5
powers <- comparePower(runSimsResults, alpha.type="fdr", alpha.nominal=0.05, 
                       stratify.by="expr", delta=0.585, filter.by="expr", 
                       strata=c(0, 5, 10, 2^(1:6) * 10, Inf), strata.filtered=1)
p1 <- summaryPower(powers)

# print some information
cat(paste0(pheno, "\n", " No Exclusion: ", round(p0[,4],2), "\n", 
           " After exclusion: ", round(p1[,4],2), "\n\n"))





                      
