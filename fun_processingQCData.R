#this fonction performe the QC of the raw data BEFORE NORMALIZATION.
processingQCData <- function(path, server){

####QC proposed by minfi
    #load library
      if (server=="distant"){
	  LIB_METH = NULL
	  }else{
	  LIB_METH = "~/share/greenwood.group/Rlibs/methylationR3.1"}

	  require(minfi, lib.loc = LIB_METH  )

      #load R script
       #load R script
      source("~/methylation450KPipeline/extractFromRGSet450k.R")
      #source("/home/data1/homeldi/greg.voisin/share/greenwood.group/Rscripts/methylation_R_scripts/extractFromRGSet450k.R")
      
    #load R objects
      load( paste(path,"RAW_DATA/RGset.RData", sep = "/"))
      pd <- pData(RGset)

    # IF sample group is empty(NA ) in the sample sheet.
      if(sum(is.na( pd$Sample_Group)) == length( pd$Sample_Group)){
         pd$Sample_Group <- rep(1, length( pd$Sample_Group))
      }

    # generate a R pdf Report from minfi
      qcReport(RGset, sampNames = pd$Sample_Name, sampGroups = pd$Sample_Group,  pdf = paste(path,"QC_REPORT","qcReport.pdf", sep = "/"))

####QC proposed by JP/Greg
    #load data and control values
    extractedData <- extractFromRGSet450k(RGset)
    #create a list for greenControl
    greenControls <- extractedData$greenControls
    #create a list for redControl
    redControls <- extractedData$redControls
    #vector of the names of kind of controls.
    controlNames <- names(greenControls)

    ##. Dye bias study
    index <- match("NEGATIVE",controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    dyebiasMatrix <- log2(greenControls.current / redControls.current)
    dyebias <- apply(dyebiasMatrix, 2, median)

    svg(paste(path,"QC_REPORT/biaisDyeBias.svg", sep= "/"))
      par(mar = c(10,5,4,1))
      boxplot(dyebiasMatrix, main= "boxplot log2(GreenSignal/RedSignal)\n for NegativeControls", las= 2, ylab= "ratio Green/Red")
    dev.off()

    ## Bisulfite conversion extraction for probe type II:
    index <-    match("BISULFITE CONVERSION II",controlNames)
    redControls.current <- redControls[[ index ]]
    bisulfite2 <- colMeans(redControls.current, na.rm = TRUE)

    svg(paste(path,"QC_REPORT/BISULFITE_CONVERSION_II.svg", sep= "/"))
      par(mar = c(10,5,4,1))
      plot(redControls.current[1,] , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(redControls.current), max(redControls.current)), main = "BISULFITE CONVERSION II",ylab = "signal", xlab = "", xaxt="n")
      axis(1, at=1:ncol(redControls.current), labels=colnames(redControls.current), padj = 1, las= 2)
        for ( i in 2: nrow(redControls.current)){
            points(redControls.current[i,], col = colors()[sample(1:657,1)], type= "l")
        }
    dev.off()

    ## Bisulfite conversion extraction for probe type I:
    index <- match("BISULFITE CONVERSION I",controlNames)
    redControls.current <- redControls[[ index ]][7:9,]
    vminRed = min(redControls.current)
    vmaxRed = max(redControls.current)
    greenControls.current <- redControls[[ index ]][1:3,]
    vminGreen = min(greenControls.current)
    vmaxGreen = max(greenControls.current)
    bisulfite1 <- colMeans(redControls.current + greenControls.current, na.rm = TRUE)

    svg(paste(path,"QC_REPORT/BISULFITE_CONVERSION_I.svg", sep= "/"))
      par(mar = c(10,5,4,1))
      plot(redControls.current[1,] , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(c(vminRed, vminGreen)), max(c(vmaxRed,vmaxGreen))), main = "BISULFITE CONVERSION I",ylab = "signal", xlab = "", xaxt="n")
      axis(1, at=1:ncol(redControls.current), labels=colnames(redControls.current), padj = 1, las= 2)
        for ( i in 2: nrow(redControls.current)){
            points(redControls.current[i,], col = colors()[sample(1:657,1)], type= "l")
        }
        for ( i in 1: nrow(greenControls.current)){
            points(greenControls.current[i,], col = colors()[sample(1:657,1)], type= "l", lty=2)
        }
    dev.off()

    ## Staining
    index <- match("STAINING",controlNames)
    sg <- greenControls[[ index ]][3, ]
    sr <- redControls[[ index ]][1, ]
    staining_value = cbind(sg,sr)

    svg(paste(path,"QC_REPORT/STAINING.svg", sep= "/"))
      par(mar = c(10,5,4,1))
      plot(staining_value[,1] , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(staining_value), max(staining_value)), main = "STAINING ",ylab = "signal", xlab = "", xaxt="n")
      axis(1, at=1:ncol(redControls.current), labels=names(sr), padj = 1, las= 2)
      points(staining_value[,2], col = colors()[sample(1:657,1)], type= "l")
    dev.off()


    ## Extension
    index <-    match("EXTENSION",controlNames)
    redControls.current     <- redControls[[index]]
    vminRed = min(redControls.current)
    vmaxRed = max(redControls.current)
    greenControls.current   <- greenControls[[index]]
    vminGreen = min(greenControls.current)
    vmaxGreen = max(greenControls.current)
    extr <- redControls.current[1:2,]
    extg <- greenControls.current[3:4,]

    svg(paste(path,"QC_REPORT/EXTENSION.svg", sep= "/"))
    par(mar = c(10,5,4,1))
    plot(redControls.current[1,] , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(c(vminRed, vminGreen)), max(c(vmaxRed,vmaxGreen))), main = "EXTENSION",ylab = "signal", xlab = "", xaxt="n")
   axis(1, at=1:ncol(redControls.current), labels=colnames(redControls.current), padj = 1, las= 2)
   for ( i in 2: nrow(redControls.current)){
       points(redControls.current[i,], col = colors()[sample(1:657,1)], type= "l")
   }
   for ( i in 1: nrow(greenControls.current)){
       points(greenControls.current[i,], col = colors()[sample(1:657,1)], type= "l", lty=2)
   }
   dev.off()


   ## Hybridization should be monitored only in the green channel
   index <- match("HYBRIDIZATION",controlNames)
   h1 <- greenControls[[index]][1, ]
   h2 <- greenControls[[index]][2, ]
   h3 <- greenControls[[index]][3, ]

   svg(paste(path,"QC_REPORT/HYBRIDIZATION.svg", sep= "/"))
      par(mar = c(10,5,4,1))
      plot(h1 , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(c(h1,h2,h3)), max(c(h1,h2,h3))), main = "HYBRIDIZATION",ylab = "signal", xlab = "", xaxt="n")
      points(h2, col = colors()[sample(1:657,1)], type= "l")
      points(h3, col = colors()[sample(1:657,1)], type= "l")
      axis(1, at=1:ncol(redControls.current), labels=colnames(redControls.current), padj = 1, las= 2)
   dev.off()


   ##Target removal should be low compared to hybridization probes
   index <- match("TARGET REMOVAL",controlNames)
   tar <- greenControls[[index]]

   svg(paste(path,"QC_REPORT/TARGET_REMOVAL.svg", sep= "/"))
     par(mar = c(10,5,4,1))
     plot(tar[1,] , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(tar), max(tar)), main = "TARGET_REMOVAL",ylab = "signal", xlab = "", xaxt="n")
     points(tar[2,], col = colors()[sample(1:657,1)], type= "l")
     axis(1, at=1:ncol(tar), labels=colnames(tar), padj = 1, las= 2)
   dev.off()

   ## Non-polymorphic probes
   index <- match("NON-POLYMORPHIC",controlNames)
   npr <- redControls[[index]][1:2, ]
   npr_means <- colMeans(npr, na.rm = TRUE)

   vminRed = min(npr)
   vmaxRed = max(npr)
   npg <- greenControls[[index]][3:4, ]
   npg_means <- colMeans(npg, na.rm = TRUE)

   vminGreen = min(npg)
   vmaxGreen = max(npg)


   svg(paste(path,"QC_REPORT/NON_POLYMORPHIC.svg", sep= "/"))
     par(mar = c(10,5,4,1))
     plot(npr[1,] , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(c(vminRed, vminGreen)), max(c(vmaxRed,vmaxGreen))), main = "NON_POLYMORPHIC",ylab = "signal", xlab = "", xaxt="n")
     axis(1, at=1:ncol(npr), labels=colnames(npr), padj = 1, las= 2)
       for ( i in 2: nrow(npr)){
           points(npr[i,], col = colors()[sample(1:657,1)], type= "l")
       }
       for ( i in 1: nrow(npg)){
           points(npg[i,], col = colors()[sample(1:657,1)], type= "l", lty=2)
       }
   dev.off()

   ## Specificity II
   index <- match("SPECIFICITY II",controlNames)
   greenControls.current <- greenControls[[index]]
   redControls.current <- redControls[[index]]
   spec2g <- colMeans(greenControls.current, na.rm = TRUE)
   spec2r <- colMeans(redControls.current, na.rm = TRUE)
   spec2ratio <- spec2g / spec2r
   spec2g <- greenControls.current
   vminGreen = min(spec2g)
   vmaxGreen = max(spec2g)
   spec2r <- redControls.current
   vminRed = min(spec2r)
   vmaxRed = max(spec2r)

   svg(paste(path,"QC_REPORT/SPECIFICITY_II.svg", sep= "/"))
     par(mar = c(10,5,4,1))
     plot(spec2g[1,] , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(c(vminRed, vminGreen)), max(c(vmaxRed,vmaxGreen))), main = "SPECIFICITY_II",ylab = "signal", xlab = "", xaxt="n")
     axis(1, at=1:ncol(spec2g), labels=colnames(spec2g), padj = 1, las= 2)
       for ( i in 2: nrow(spec2g)){
           points(spec2g[i,], col = colors()[sample(1:657,1)], type= "l")
       }
       for ( i in 1: nrow(spec2r)){
           points(spec2r[i,], col = colors()[sample(1:657,1)], type= "l", lty=2)
       }
   dev.off()

    ## Specificity I
    index <- match("SPECIFICITY I", controlNames)
    greenControls.current <- greenControls[[index]][1:3,]
    redControls.current <- redControls[[index]][7:9,]
    spec1g <- greenControls.current
    spec1r <- redControls.current
    greenControls.current <- greenControls[[index]][1:3,]
    redControls.current <- redControls[[index]][1:3,]
    ratio1 <- colMeans(redControls.current, na.rm = TRUE) /
    colMeans(greenControls.current, na.rm = TRUE)
    greenControls.current <- greenControls[[index]][7:9,]
    redControls.current <- redControls[[index]][7:9,]
    ratio2 <- colMeans(greenControls.current, na.rm = TRUE) /
    colMeans(redControls.current, na.rm = TRUE)
    spec1ratio <- (ratio1 + ratio2) / 2

    svg(paste(path,"QC_REPORT/SPECIFICITY_I.svg", sep= "/"))
      par(mar = c(10,5,4,1))
      plot(spec1ratio , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(spec1ratio), max(spec1ratio)), main = "SPECIFICITY_I",ylab = "signal", xlab = "", xaxt="n")
      axis(1, at=1:length(spec1ratio), labels=colnames(spec1ratio), padj = 1, las= 2)
    dev.off()

    ## Normalization probes:
    index <- match(c("NORM_A"), controlNames)
    normA <- colMeans(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_T"), controlNames)
    normT <- colMeans(redControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_C"), controlNames)
    normC <- colMeans(greenControls[[index]], na.rm = TRUE)
    index <- match(c("NORM_G"), controlNames)
    normG <- colMeans(greenControls[[index]], na.rm = TRUE)
    dyebias2 <- (normC + normG)/(normA+normT)

     svg(paste(path,"QC_REPORT/NORMALIZATION_PROBES.svg", sep= "/"))
       par(mar = c(10,5,4,1))
       plot(dyebias2 , col = colors()[sample(1:657,1)], type= "l" , ylim = c(min(dyebias2), max(dyebias2)), main = "NORMALIZATION_PROBES",ylab = "signal", xlab = "", xaxt="n")
       axis(1, at=1:length(dyebias2), labels=colnames(dyebias2), padj = 1, las= 2)
     dev.off()

	## compute the MAD score
	load(paste(path,"RAW_DATA/raw_Bvalue.RData", sep = "/"))
	madscore = apply(raw_Bvalue,2, mad)
	madscore_ordered = madscore[order(madscore)]
	svg(paste(path,"QC_REPORT/MADSCORE.svg", sep= "/"))
       par(mar = c(10,5,4,1))
       plot(madscore  ,col = "grey",cex= 0.4, ylim = c(min(madscore), max(madscore)), main = "MADSCORE",ylab = "MADScore", xlab = "", xaxt="n")
       points(match(names(madscore_ordered)[1:10],names(madscore)),madscore_ordered[1:10], col= "red")
	   text(match(names(madscore_ordered)[1:10],names(madscore)),madscore_ordered[1:10], labels = names(madscore_ordered)[1:10], cex= 0.6, col = "red")
	   axis(1, at=1:length(madscore), labels=colnames(madscore), padj = 1, las= 2)
	dev.off()


	## compute the MAD score
	load(paste(path,"RAW_DATA/detection_pvalue.RData", sep = "/"))
	nb0.01 =  apply(detection_pvalue, 2, function(x) sum(x < 0.01))
	nb0.05 =  apply(detection_pvalue, 2, function(x) sum(x < 0.05))
	detected_sum = nb0.01+nb0.05
	detected_sum_ordered = detected_sum[order(detected_sum)]

	valueX = nb0.01[match(names(detected_sum_ordered)[1:10],names(nb0.01))]
	valueY = nb0.05[match(names(detected_sum_ordered)[1:10],names(nb0.05))]



	svg(paste(path,"QC_REPORT/Detected_CpG.svg", sep= "/"))
       par(mar = c(10,5,4,1))
       plot(nb0.01, nb0.05 ,col = "grey",cex= 0.4, ylim = c(min(nb0.01), max(nb0.01)), main = "Detected Pvalue",ylab = "detected(0.05)", xlab = "detected(0.01)")
       points(valueX,valueY, col= "red")
	   text(valueX,valueY, labels = names(valueX)[1:10], cex= 0.6, col = "red",pos= 4)
     dev.off()



    #create a file to summarize the QC
    #compute some means
 extr_mean    = colMeans(extr, na.rm = TRUE)
 extg_mean    = colMeans(extg, na.rm = TRUE)
 spec2r_mean  = colMeans(spec2r, na.rm = TRUE)
 spec2g_mean  = colMeans(spec2g, na.rm = TRUE)
 spec1r_mean  = colMeans(spec1r, na.rm = TRUE)
 spec1g_mean  = colMeans(spec1g, na.rm = TRUE)
 tar_means    = colMeans(tar, na.rm = TRUE)

 #compile data
  QC_value_PARAMETER <- cbind(dyebias, bisulfite2, bisulfite1,staining_value, extr_mean, extg_mean, h1,h2,h3, tar_means,npr_means, npg_means, spec2r_mean, spec2g_mean, spec2ratio , spec1r_mean, spec1g_mean, spec1ratio, dyebias2, madscore, nb0.01, nb0.05)

  QC_value_PARAMETER = apply(QC_value_PARAMETER, 2, function(x) {return(round(x,2))})

 colnames(QC_value_PARAMETER) <- c("dyeBiasMedian" , "bisulfite2", "bisulfite1","staining_green", "staining_red", "extension_green", "extension_red", "hybrid_1", "hybrid_2", "hybrid_3","target_removal", "non_polymorp_red","non_polymorp_green", "specificII_red","specificII_green", "specificII_ratio", "specificI_red","specificI_green","specificI_ratio","dyebias2", "MADScore", "Detected(0.01)","Detected(0.05)")



  #compute the Statut for 6 parameters expected to be > 10000 : bisulfite 1, bisulfite 2, specific 1,2 green and red.
  bisulfite2_STATUT                     <- rep(6, length(bisulfite2))
  bisulfite2_STATUT[bisulfite2 < 10000] <- 5
  bisulfite2_STATUT[bisulfite2 < 5000]  <- 4
  bisulfite2_STATUT[bisulfite2 < 4000]  <- 3
  bisulfite2_STATUT[bisulfite2 < 3000]  <- 2
  bisulfite2_STATUT[bisulfite2 < 2000]  <- 1
  bisulfite2_STATUT[bisulfite2 < 1000]  <- 0

  bisulfite1_STATUT                     <- rep(6, length(bisulfite1))
  bisulfite1_STATUT[bisulfite1 < 10000] <- 5
  bisulfite1_STATUT[bisulfite1 < 5000]  <- 4
  bisulfite1_STATUT[bisulfite1 < 4000]  <- 3
  bisulfite1_STATUT[bisulfite1 < 3000]  <- 2
  bisulfite1_STATUT[bisulfite1 < 2000]  <- 1
  bisulfite1_STATUT[bisulfite1 < 1000]  <- 0

  spec2r_STATUT                      <- rep(6, length(spec2r_mean))
  spec2r_STATUT[spec2r_mean < 10000] <- 5
  spec2r_STATUT[spec2r_mean < 5000]  <- 4
  spec2r_STATUT[spec2r_mean < 4000]  <- 3
  spec2r_STATUT[spec2r_mean < 3000]  <- 2
  spec2r_STATUT[spec2r_mean < 2000]  <- 1
  spec2r_STATUT[spec2r_mean < 1000]  <- 0

  spec2g_STATUT                      <- rep(6, length(spec2g_mean))
  spec2g_STATUT[spec2g_mean < 10000] <- 5
  spec2g_STATUT[spec2g_mean < 5000]  <- 4
  spec2g_STATUT[spec2g_mean < 4000]  <- 3
  spec2g_STATUT[spec2g_mean < 3000]  <- 2
  spec2g_STATUT[spec2g_mean < 2000]  <- 1
  spec2g_STATUT[spec2g_mean < 1000]  <- 0

  spec1r_STATUT                      <- rep(6, length(spec1r_mean))
  spec1r_STATUT[spec1r_mean < 10000] <- 5
  spec1r_STATUT[spec1r_mean < 5000]  <- 4
  spec1r_STATUT[spec1r_mean < 4000]  <- 3
  spec1r_STATUT[spec1r_mean < 3000]  <- 2
  spec1r_STATUT[spec1r_mean < 2000]  <- 1
  spec1r_STATUT[spec1r_mean < 1000]  <- 0

  spec1g_STATUT                     <- rep(6, length(spec1g_mean))
  spec1g_STATUT[spec1g_mean < 10000] <- 5
  spec1g_STATUT[spec1g_mean < 5000]  <- 4
  spec1g_STATUT[spec1g_mean < 4000]  <- 3
  spec1g_STATUT[spec1g_mean < 3000]  <- 2
  spec1g_STATUT[spec1g_mean < 2000]  <- 1
  spec1g_STATUT[spec1g_mean < 1000]  <- 0


  rank_MADScore = 1:length(madscore_ordered)
  names(rank_MADScore) = names(madscore_ordered)
  rank_MADScore = rank_MADScore[names(bisulfite2)] # use names(bisulfite2) to order the MASCore in the same order (  the order function don't give a expected rank)

  rank_detected_sum = 1:length(detected_sum_ordered)
  names(rank_detected_sum) = names(detected_sum_ordered)
  rank_detected_sum = rank_detected_sum[names(bisulfite2)]

  QC_value_STATUT <-cbind(bisulfite2_STATUT,bisulfite1_STATUT, spec2r_STATUT,spec2g_STATUT,spec1r_STATUT,spec1g_STATUT )
  STATUT_VALUE = apply(QC_value_STATUT,1,sum)
  QC_value_STATUT =cbind(QC_value_STATUT, rank_MADScore, rank_detected_sum)


	svg(paste(path,"QC_REPORT/STATUT_VALUE.svg", sep= "/"))
       par(mar = c(10,5,4,1))
       plot(STATUT_VALUE ,col = "grey",cex= 0.4, ylim = c(min(STATUT_VALUE), max(STATUT_VALUE)), main = "SCORE_QC",ylab = "SCORE_QC", xlab = "", xaxt="n")
       axis(1, at=1:length(STATUT_VALUE), labels=colnames(STATUT_VALUE), padj = 1, las= 2)
     dev.off()



  #combine data.

  QC_value_PARAMETER = cbind(sample_ID = names(dyebias),STATUT_VALUE,QC_value_PARAMETER,QC_value_STATUT)

 #   #read the covariables files to associate Covariable to the result of the QC.
#     if (file.exists(paste(path,"COVARIABLE_DATA/CovariableModified.csv", sep= "/"))){
#       covariables_file <- read.csv2(paste(path,"COVARIABLE_DATA/CovariableModified.csv", sep= "/"), sep= "," ,header = T)
#     }else{
#       covariables_file <- read.csv2(paste(path,"COVARIABLE_DATA/Covariable.csv", sep= "/"), sep= "," ,header = T)
#     }#end of if
# 	Sample_ID = paste(covariables_file$Sentrix_ID,covariables_file$Sentrix_Position, sep= "_")
# 	covariables_file = cbind(Sample_ID, covariables_file)
# 
#   #for test
#   #save(QC_value_PARAMETER, file = "QC_value_PARAMETER.RData")
#   #save(covariables_file,   file = "covariables_file.RData")
# 
#    QC_value_PARAMETER = merge(QC_value_PARAMETER, covariables_file , by.x =  "sample_ID", by.y = "Sample_ID")
#    QC_value_PARAMETER = QC_value_PARAMETER[order(as.numeric(as.vector(QC_value_PARAMETER$STATUT_VALUE)) ,decreasing = F),]

   # write the result of the QC in order of the QC value ( more the QC value is hight, more the QC is better
   write.csv(QC_value_PARAMETER, paste(path, "QC_REPORT/QC_value_summary.csv", sep = "/"),quote = F,  eol = "\n", row.names = F)


  # process the html report.
  library(rmarkdown)
  # exige blade 7 becauce PanDoc
  setwd(paste(path,"QC_REPORT", sep= "/"))
  render("QCReport.Rmd", "html_document")

}#end of processingQCData
