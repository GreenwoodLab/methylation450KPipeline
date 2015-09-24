source(paste(Sys.getenv("PIPELINE_450K"), "/fun_initialize_analysis.R", sep="")) 

print("..Initialize the analysis: create the folder")

    #get the arguments from batch command line
    args <- commandArgs(trailingOnly = TRUE)

    #process the process of IdatData before inspection of QC
    initialize_analysis(args[1])

print("..End of the creation of the Folder")
