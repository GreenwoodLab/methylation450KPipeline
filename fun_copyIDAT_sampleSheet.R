copyIDAT_sampleSheet <- function(path1,path0){
  path1 =  paste(path1, "IDAT", sep = "/")
  system(paste("cp",paste(path0,"*.idat", sep = "/"), path1 , sep =" "))
  system(paste("cp", paste(path0,"*.csv", sep = "/"), path1 , sep = " "))
}
