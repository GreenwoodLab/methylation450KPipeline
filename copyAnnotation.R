#script for methylation pipeline analysis
#february 2014
#voisin greg - Greenwood lab- ldi

################################################################################
source("fun_copyAnnotation.R")

print(paste("Beginning of copy of Annotation :\n", Sys.time(),sep= ""))

    #get the arguments from batch command line
    args <- commandArgs(trailingOnly = TRUE)

    #process copy of annotation data
    copyAnnotation(args[1])

print("..The annotation file  are copied in ANNOTATION folder of the analysis")
