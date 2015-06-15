
print(paste("Beginning of the analysis :\n", Sys.time(),sep= ""))
 
print("..Initialize the analysis: create the folder")

    #get the arguments from batch command line
    args <- commandArgs(trailingOnly = TRUE)

    #process the process of IdatData before inspection of QC
    initialize_analysis(args[1])

print("..End of the creation of the Folder")
