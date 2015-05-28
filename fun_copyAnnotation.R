#Copy all annotation data
#  input: path1: folder analysis
copyAnnotation <- function(path1, path0){
  path0 = "~/ANNOTATION_FOR450K_PIPELINE"
  path1 = paste(path1, "ANNOTATION_DATA", sep = "/")
  system(paste("cp",paste(path0,"*.*", sep = "/"), path1 , sep = " "))
}
