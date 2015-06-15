#script for methylation pipeline analysis
#february 2014
#voisin greg - Greenwood lab- ldi

################################################################################
print(paste("Beginning of copy of IDAT file and samplesheet :\n", Sys.time(),sep= ""))

    #get the arguments from batch command line
    args <- commandArgs(trailingOnly = TRUE)

    #process the copy of Idat files 
    copyIDAT_sampleSheet(args[1], args[2])

print("..The idat files and csv are copied in IDAT folder of the analysis")
