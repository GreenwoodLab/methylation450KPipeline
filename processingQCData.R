source("fun_processingQCData.R")

print(paste("Beginning the Quality Control of the Idat data :\n", Sys.time(),sep= ""))
    
    #get the arguments from batch command line
    args <- commandArgs(trailingOnly = TRUE)

    #process the process of IdatData before inspection of QC
    processingQCData(args[1],args[2])
    
print("..The QC of Idat files have been processed")
