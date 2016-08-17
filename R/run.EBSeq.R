#
# A. Gaye 16th august 2016 #
# A new function to use with the library POWER* to enable power analysis using 'EBSeg'** 
# a method more suited for isoform level differential expresison analyis.
# *  Bioinformatics. 2015 Jan 15;31(2):233-41. doi: 10.1093/bioinformatics/btu640
# ** Bioinformatics (2013) 29 (8): 1035-1043. doi: 10.1093/bioinformatics/btt087 
#

run.EBSeq <- function(dat){
  
  # generate the matrix of counts, the vectors of gene names and isoform names and store all 3 objects in a list
  IsoMat <- dat$counts
  colnames(IsoMat) <- NULL
  rownames(IsoMat) <- gsub(" ", ".", paste0(dat$iso$gene_id, "|", dat$iso$transcript_id), fixed=TRUE)
  # generate the matrix object
  IsoNames <- gsub(" ", ".", paste0(dat$iso$gene_id, "|", dat$iso$transcript_id), fixed=TRUE)
  IsosGeneNames <- gsub(" ", ".", dat$iso$gene_id, fixed=TRUE)
  IsoList <- list("IsoMat"=IsoMat, "IsoNames"=IsoNames, "IsosGeneNames"=IsosGeneNames)
  
  # get normalization factors using the TMM method instead of the default EBSeq's median normalization 
  #cat("normalize\n")
  IsoSizes <- calcNormFactors(IsoMat)
  
  # prior parameters for different uncertainty groups
  NgList <- GetNg(IsoNames, IsosGeneNames, 3)
  IsoNgTrun <- NgList$IsoformNgTrun
  
  # identify cases and controls
  status <- dat$designs
  status[which(status==1)] <- "C1"
  status[which(status==0)] <- "C2"
  
  # run DE analysis and extract results
  filter1 <- 1
  filter2 <- 0
  numiter <- 2
  #cat("test\n")
  IsoEBOut <- EBTest(Data=IsoMat, NgVector=IsoNgTrun, 
                     Conditions=as.factor(status), 
                     sizeFactors=IsoSizes, maxround=numiter, Qtrm=filter1, QtrmCut=filter2)
  
  # check that hyper-parameter estimations are converged if not increase number of iterations by 5
  # Do this for up to 50 iterations if convergence is still not achieve print a warning and carry on
  for(i in 1:18){
    # check if convergence achieved (i.e. if the 2 last estimates are equal) break out of loop
    lstEst <- length(IsoEBOut$Alpha)
    forlstEst <- length(IsoEBOut$Alpha)-1
    if(round(IsoEBOut$Alpha[lstEst], 6) == round(IsoEBOut$Alpha[lstEst], 6) & 
         round(IsoEBOut$Beta[lstEst], 6) == round(IsoEBOut$Beta[lstEst], 6)){
      message(paste0("hyper-parameter estimates converged after ", numiter, " iterations!"))
      break
    }else{
      numiter <- numiter+1
      IsoEBOut <- EBTest(Data=IsoMat, NgVector=IsoNgTrun, 
                         Conditions=as.factor(status), 
                         sizeFactors=IsoSizes, maxround=numiter, Qtrm = filter1, QtrmCut = filter2)
      if(numiter == 20) {message("PLEASE VERIFY: HYPER-PARAMETER ESTIMATES DID NOT CONVERGENCE AFTER 20 ITERATIONS!")}
    }
  }

  # get the posterior probabilities of being DE or EE (equally expressed)
  IsoEBDERes <- GetPPMat(IsoEBOut)
  
  # extract the posterior fold changes and append them to the PP matrix above
  res <- PostFC(IsoEBOut)
  PostLogFC <- log2(res$PostFC)
  isoform <- rownames(IsoEBDERes)
  IsoEBDERes <- cbind(isoform, PostLogFC, IsoEBDERes)
  
  # generate final results to return
  fdr <- (1-(as.numeric(IsoEBDERes[,"PPDE"])))
  fdrOdered <- sort(fdr)
  ix <- sort(fdr, index.return = TRUE)$ix
  result <- data.frame(geneIndex = ix, pval = fdrOdered*length(fdr), fdr = fdrOdered)
  res <- result[order(result$geneIndex), ]
  return(res)
  
}

