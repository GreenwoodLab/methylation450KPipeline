adjusted.Bvalue <- function(matrixValue){
                    for ( i in 1:ncol(matrixValue)){ matrixValue[,i] <- (matrixValue[,i]*(1000-1)+0.5)/1000
                     return(matrixValue)
                     }
                   }#end of adjusted.Bvalue
