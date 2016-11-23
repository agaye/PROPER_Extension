#
# A. Gaye 16th August 2016 #
# Example on how to use the new and amended functions 
# to run power analysis for isoform level expression data.
# Please consider looking into the documentation of the 
# libraries PROPER and EBSeq for details about the arguments 
# the functions called in this script.
#

# load required R libraries 
library(PROPER)
library(EBSeq)
library(edgeR)
source("run.EBSeq.R")
source("runSims_ag.R")
source("plotPower_ag.R")
source("plotFDR_ag.R")

# load data 
load("toy_data.RData")
load("AveLogCPM.RData")
load("logDispersion.RData")
load("libCount.RData")
load("logFC.RData")


# proportion of DE genes
propDE <- 0.05

# set the simulations
mm <- RNAseq.SimOptions.2grp(ngenes=dim(toy_data)[1], lBaselineExpr=AveLogCPM, 
                             lOD=logDispersion, seqDepth=libCount, 
                             lfc=logFC, p.DE=propDE, sim.seed=123)


# run simulation and test for DE
runSimsResults <- runSims_ag(Nreps=4, sim.opts=mm, DEmethod="EBSeq", isoforms=toy_data[,c(1,2)], nsims=20)

# evakuate power across all simulation runs
powers <- comparePower(runSimsResults, alpha.type="fdr", alpha.nominal=0.05, 
                       stratify.by="expr", filter.by="none", 
                       strata=c(0, 5, 10, 2^(1:6) * 10, Inf), strata.filtered=0)

# plot results of the power analysis
pdf("myPlots.pdf", width=10, height=7)
par(mfcol=c(1,2), oma = c(0, 0, 2, 0))
plotPower_ag(powers, error.bar=FALSE)
plotFDR(powers, error.bar=FALSE)
dev.off()


