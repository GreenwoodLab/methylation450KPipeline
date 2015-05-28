## Extraction of the Control matrix
buildControlMatrix450k <- function(extractedData) {

    greenControls <- extractedData$greenControls
    redControls <- extractedData$redControls
    controlNames <- names(greenControls)


    ## Dye bias:
    index <- match("NEGATIVE",controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    dyebiasMatrix <- log2(greenControls.current / redControls.current)
    dyebias <- apply(dyebiasMatrix, 2, median)

    ## Bisulfite conversion extraction for probe type II:
    index <-    match("BISULFITE CONVERSION II",controlNames)
    redControls.current <- redControls[[ index ]]
    bisulfite2 <- colMeans(redControls.current, na.rm = TRUE)

    ## Bisulfite conversion extraction for probe type I:
    index <- match("BISULFITE CONVERSION I",controlNames)
    redControls.current <- redControls[[ index ]][7:9,]
    greenControls.current <- redControls[[ index ]][1:3,]
    bisulfite1 <- colMeans(redControls.current + greenControls.current, na.rm = TRUE)

    ## Staining
    index <- match("STAINING",controlNames)
    sg <- greenControls[[ index ]][3, ]
    sr <- redControls[[ index ]][1, ]

    ## Extension
    index <-    match("EXTENSION",controlNames)
    redControls.current     <- redControls[[index]]
    greenControls.current   <- greenControls[[index]]
    extr <- redControls.current[1:2,]
    extg <- greenControls.current[3:4,]

    ## Hybridization should be monitored only in the green channel
    index <- match("HYBRIDIZATION",controlNames)
    h1 <- greenControls[[index]][1, ]
    h2 <- greenControls[[index]][2, ]
    h3 <- greenControls[[index]][3, ]

    ## Target removal should be low compared to hybridization probes
    index <- match("TARGET REMOVAL",controlNames)
    tar <- greenControls[[index]]

    ## Non-polymorphic probes
    index <- match("NON-POLYMORPHIC",controlNames)
    npr <- redControls[[index]][1:2, ]
    npg <- greenControls[[index]][3:4, ]

    ## Specificity II
    index <- match("SPECIFICITY II",controlNames)
    greenControls.current <- greenControls[[index]]
    redControls.current <- redControls[[index]]
    spec2g <- colMeans(greenControls.current, na.rm = TRUE)
    spec2r <- colMeans(redControls.current, na.rm = TRUE)
    spec2ratio <- spec2g / spec2r
    spec2g <- greenControls.current
    spec2r <- redControls.current

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

    model.matrix <- cbind(
        bisulfite1, bisulfite2, t(extg), t(extr), h1, h2,h3, sg, sr, t(npg),
        t(npr), t(tar), t(spec1g), t(spec1r), t(spec2g), t(spec2r),ratio1,
        spec1ratio, spec2ratio, ratio2, normA, normC, normT, normG, dyebias2)

    oobG <- extractedData$oob$greenOOB
    oobR <- extractedData$oob$redOOB
    model.matrix <- cbind(model.matrix, t(oobG[c(1,50,99),]),oobG[50,]/oobR[50,])

    ## Imputation
    for (colindex in 1:ncol(model.matrix)) {
        column <- model.matrix[,colindex]
        column[is.na(column)]    <- mean(column, na.rm = TRUE)
        model.matrix[ , colindex] <- column
    }

    ## Scaling
    model.matrix <- scale(model.matrix)

    ## Fixing outliers
    model.matrix[model.matrix > 3] <- 3
    model.matrix[model.matrix < (-3)] <- -3

    ## Rescaling
    model.matrix <- scale(model.matrix)

    return(model.matrix)
}
