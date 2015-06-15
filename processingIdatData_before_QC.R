source("methylation450KPipeline/fun_processingIdatData_before_QC.R")

print(paste("Beginning the Idat data process of :", Sys.time(),sep= ""))
    
    #get the arguments from batch command line
    args <- commandArgs(trailingOnly = TRUE)

    #process the process of IdatData before inspection of QC
    processingIdatData_before_QC(args[1], args[2])
    
print("..The Idat file have been processed")
