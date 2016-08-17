#
# A. Gaye 16th august 2016 #
# This is an amended version of the function 'runSim' from the library POWER* to enable 
# power analysis using 'EBSeg'** a method more suited for isoform level differential expresison analyis.
# *  Bioinformatics. 2015 Jan 15;31(2):233-41. doi: 10.1093/bioinformatics/btu640
# ** Bioinformatics (2013) 29 (8): 1035-1043. doi: 10.1093/bioinformatics/btt087 
#

runSims_ag <- function (Nreps = c(3, 5, 7, 10), nsims = 100, sim.opts, 
                        DEmethod = c("EBSeq","edgeR", "DSS", "DESeq"), isoforms = NULL, verbose = TRUE) 
{
  DEmethod = match.arg(DEmethod)
  n1 = n2 = max(Nreps)
  set.seed(sim.opts$sim.seed)
  pvalue = fdrs = xbar = array(NA, dim = c(sim.opts$ngenes, 
                                           length(Nreps), nsims))
  DEids = lfcs = NULL
  for (i in 1:nsims) {
    if (verbose) 
      cat("Simulation number", i, "\n")
    sim.opts$sim.seed = sim.opts$sim.seed + 1
    sim.opts = PROPER:::update.RNAseq.SimOptions.2grp(sim.opts)
    dat.sim.big = PROPER:::simRNAseq(sim.opts, n1, n2)
    DEids[[i]] = dat.sim.big$DEid
    lfcs[[i]] = dat.sim.big$simOptions$lfc
    for (j in seq(along = Nreps)) {
      Nrep = Nreps[j]
      idx = c(1:Nrep, max(Nreps) + (1:Nrep))
      this.design = dat.sim.big$designs[idx]
      this.X = dat.sim.big$counts[, idx]
      this.X[which(is.na(this.X))] <- 0
      this.simOpts = sim.opts
      ss = rowSums(this.X)
      ix.valid = ss > 0
      this.X.valid = this.X[ix.valid, , drop = FALSE]
      iso.valid = isoforms[ix.valid, , drop = FALSE]
      if(DEmethod == "EBSeq"){
        if(is.null(isoforms)){stop("Please provide a 2-column table with the names of the isoforms!")}
        data0 = list(counts = this.X.valid, designs = this.design, iso=iso.valid)
      }else{
        data0 = list(counts = this.X.valid, designs = this.design)
      }
      if (DEmethod == "EBSeq") 
        print(dim(data0$counts))
        res1 = run.EBSeq(data0)
      if (DEmethod == "edgeR") 
        res1 = run.edgeR(data0)
      if (DEmethod == "DESeq") 
        res1 = run.DESeq(data0)
      if (DEmethod == "DSS") 
        res1 = run.DSS(data0)
      pval = fdr = rep(1, nrow(this.X))
      X.bar1 = rep(0, nrow(this.X))
      pval[ix.valid] = res1[, "pval"]
      fdr[ix.valid] = res1[, "fdr"]
      sizeF = colSums(data0$count)
      sizeF = sizeF/median(sizeF)
      X.bar1[ix.valid] = rowMeans(sweep(data0$count, 2, 
                                        sizeF, FUN = "/"))
      pvalue[, j, i] = pval
      fdrs[, j, i] = fdr
      xbar[, j, i] = X.bar1
    }
  }
  list(pvalue = pvalue, fdrs = fdrs, xbar = xbar, DEid = DEids, 
       lfcs = lfcs, Nreps = Nreps, sim.opts = sim.opts)
}
