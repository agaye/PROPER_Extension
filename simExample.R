#
# A. Gaye 16th August 2016 #
# Example on how to use the new and amended functions 
# to run power analysis for isoform level expression data.
# Please look into the R library 'PROPER'* for the arguments
# of the main function 'RNAseq.SimOptions.2grp'
# * Bioinformatics. 2015 Jan 15;31(2):233-41. doi: 10.1093/bioinformatics/btu640
# N.B: This is just an example and the results are just an illustration and 
# have no real significance
#

setwd("/home/gayea/Publications/gExpression01/DEpower/PROPER4EBSeq/POWERextension")
# load required R libraries
library(PROPER)
library(edgeR)
library(EBSeq)

# load the expresison and phenotype data
load("exprDat.RData")
load("phenoDat.RData")
t <- exprDat
r0 <-  phenoDat

# filter out genes that did not achieve 1 count per million in at least 3 samples
y <- DGEList(counts=t[,3:dim(t)[2]], genes=t[, 1])
keep0 <- rowSums(cpm(y)>1) >= 3
y <- y[keep0,]
isonames <- t[keep0,c(1,2)]

# identify the cases and choose randomly an equal number of cases and controls
# This is because in the current version of POWER to two groups must be of equal size.
cases <- as.character(r0[which(r0[,"OUTCOME"] == 1), 1])
cases <- cases[which(cases %in% colnames(y$counts))]
controls <- as.character(r0[which(r0[,"OUTCOME"] == 0), 1])
controls <- controls[which(controls %in% colnames(y$counts))]
tokeep <- c(cases, sample(controls, length(cases)))
y <- y[,tokeep]
r0 <- r0[which(r0$ID %in% tokeep),]

# The model to fit
y <- calcNormFactors(y)
lp <- as.formula("~OUTCOME")
design <- model.matrix(lp, data=r0)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)

# log average CPM
aveCPM <- y$AveLogCPM

# log tagwise dispersion
logDisp <- log(y$tagwise.dispersion)

# sequence depth (total count for each library i.e. of each of the samples)
seqD <- y$samples$lib.size

# log fold change of DE genes : we used the log fold changes of genes
# genes differentially expressed in our example analysis.
# see 'PROPER' documentation for other options
logFC <- c(-3.00,-3.86,-1.18,-2.51,1.52,-0.74,-3.40,-2.22)

# proportion of DE genes
propDE <- 0.01

# generate the object that hold the parameters of the data
mm <- RNAseq.SimOptions.2grp(ngenes=dim(y$counts)[1], lBaselineExpr=aveCPM, 
                             seqDepth=seqD, lOD=logDisp, 
                             p.DE=propDE, lfc=logFC, sim.seed=123)

# run the simulation
numreps <- length(cases)
source("run.EBSeq.R")
source("runSims_ag.R")
runSimsResults <- runSims_ag(Nreps=c(numreps), sim.opts=mm, DEmethod="EBSeq", isoforms=isonames, nsims=50)
save(runSimsResults, file="runSimsResults.Rdata")
