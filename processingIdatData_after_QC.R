processingIdatData_after_QC <- function(path, server, cutoff=23){
 #load library
 if (server=="distant"){
	  LIB_METH = NULL
	  }else{
	  LIB_METH = "~/share/greenwood.group/Rlibs/methylationR3.1"}
	  
  #identify failed samples
    tmp                   <-read.csv(file.path(path,"QC_REPORT/QC_value_summary.csv"),as.is=T)
    exclusionSampleVector <- tmp$sample_ID[which(tmp$STATUT_VALUE<cutoff)]
    nb_exclusion          <-length(exclusionSampleVector)
    
   require(minfi, lib.loc = LIB_METH )

  if(nb_exclusion ==0){
  sink(file.path(path,"DOCUMENTATION/readme.txt"),append = TRUE)
  	cat("The threshold for the current analysis was set at", cutoff, "\n")
  	cat("All samples were above threshold therefore no samples were excluded from current analysis /n")
  sink()
    stop("no sample are excluded of the analysis")
  }else{
    targets <- read.450k.sheet(paste(path,"IDAT", sep= "/"))
    sink(file.path(path,"DOCUMENTATION/readme.txt"),append = TRUE)
    	 cat("The threshold for the current analysis was set at", cutoff, "\n")
  	 cat("The following samples had a statut value below the cutoff and were excluded from current analysis:\n")
  	 print(exclusionSampleVector)
    sink()
    
    exclusion_pos = NULL
    for (i in 1:nb_exclusion){
      exclusion_pos  = append(exclusion_pos, grep(exclusionSampleVector[i], targets$Basename))
    }
    targets_after_QC <- targets[-exclusion_pos,]
    RGset_after_QC   <- read.450k.exp(targets = targets_after_QC)

    #write the Fun_Norm data in normalize data folder
  FunNorm_data_bckdyeAdj_after_QC       <- preprocessFunnorm(RGset_after_QC,sex = NULL,, verbose=T)           
  FunNorm_data_NobckdyeAdj_after_QC     <- preprocessFunnorm(RGset_after_QC,sex = NULL, verbose=T, bgCorr = FALSE, dyeCorr = FALSE)          
  FunNorm_Bvalue_bckdyeAdj_after_QC     <- getBeta(FunNorm_data_bckdyeAdj_after_QC)                     
  FunNorm_Bvalue_NobckdyeAdj_after_QC   <- getBeta(FunNorm_data_NobckdyeAdj_after_QC)                         
   save(FunNorm_Bvalue_bckdyeAdj_after_QC, file = paste(path,"FUN_NORM_DATA/FunNorm_Bvalue_bckdyeAdj_after_QC.RData", sep = "/"))
   save(FunNorm_Bvalue_NobckdyeAdj_after_QC, file = paste(path,"FUN_NORM_DATA/FunNorm_Bvalue_after_QC.RData", sep = "/"))
   
    }
}
