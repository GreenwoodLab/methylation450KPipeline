extractFromRGSet450k <- function(rgSet) {
    rgSet <- updateObject(rgSet)
    controlType <- c("BISULFITE CONVERSION I",
                     "BISULFITE CONVERSION II",
                     "EXTENSION",
                     "HYBRIDIZATION",
                     "NEGATIVE",
                     "NON-POLYMORPHIC",
                     "NORM_A",
                     "NORM_C",
                     "NORM_G",
                     "NORM_T",
                     "SPECIFICITY I",
                     "SPECIFICITY II",
                     "TARGET REMOVAL",
                     "STAINING")
    MSet.raw <- preprocessRaw(rgSet)
    r <- getRed(rgSet)
    g <- getGreen(rgSet)
    meth <- getMeth(MSet.raw)
    unmeth <- getUnmeth(MSet.raw)
    beta <- getBeta(MSet.raw)
    m <- getM(MSet.raw)
    cn <- meth + unmeth

    ## Extraction of the controls
    greenControls = vector("list", length(controlType))
    redControls = vector("list", length(controlType))
    names(greenControls) <- controlType
    names(redControls) <- controlType

    for (i in 1:length(controlType)) {
        if (controlType[i] != "STAINING") {
            ctrlAddress <- getControlAddress(
                rgSet, controlType = controlType[i])
        } else {
            ctrlAddress <- getControlAddress(
                rgSet, controlType = controlType[i])[c(2,3,4,6)]
        }
        redControls[[i]] = r[ctrlAddress,]
        greenControls[[i]] = g[ctrlAddress,]

    }

    ## Extraction of the undefined negative control probes
    locusNames  <- getManifestInfo(rgSet, "locusNames")
    TypeI.Red   <- getProbeInfo(rgSet, type = "I-Red")
    TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")

    numberQuantiles <- 100
    probs <- 1:numberQuantiles/100

    ## Function to extract quantiles
    extractQuantiles <- function(matrix, probs) {
        result <- apply(matrix, 2, function(x) quantile(x, probs = probs, na.rm = TRUE))
        result
    } ## KDH: do we have a function in matrixStats?

    greenOOB <- rbind(getGreen(rgSet)[TypeI.Red$AddressA,], getGreen(rgSet)[TypeI.Red$AddressB,])
    redOOB   <- rbind(getRed(rgSet)[TypeI.Green$AddressA,], getRed(rgSet)[TypeI.Green$AddressB,])

    greenOOB <- extractQuantiles(greenOOB, probs = probs)
    redOOB   <- extractQuantiles(redOOB,   probs = probs)
    oob      <- list(greenOOB = greenOOB,  redOOB = redOOB)

    ## Defining the Type I, II Green and II Red probes:
    probesI  <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "I")
    probesII <- getProbeInfo(IlluminaHumanMethylation450kmanifest, type = "II")


    ## Chr probes:

    locations   <- getLocations(rgSet)
    chrs        <- as.character(seqnames(locations))
    names(chrs) <- names(locations)
    autosomal   <- names(chrs[chrs %in% paste0("chr", 1:22)])
    chrY        <- names(chrs[chrs %in% paste0("chr", "Y")])
    chrX        <- names(chrs[chrs %in% paste0("chr", "X")])
    probesIGrn   <- intersect( probesI$Name[probesI$Color == "Grn"], autosomal )
    probesIRed   <- intersect( probesI$Name[probesI$Color == "Red"], autosomal )
    probesII     <- intersect( probesII$Name, autosomal )
    uProbeNames <- rownames(beta)
    uProbesIGrn <- intersect(uProbeNames, probesIGrn)
    uProbesIRed <- intersect(uProbeNames, probesIRed)
    uProbesII   <- intersect(uProbeNames, probesII)
    uProbesX    <- intersect(uProbeNames, chrX)
    uProbesY    <- intersect(uProbeNames, chrY)
    indicesIGrn <- match(uProbesIGrn, uProbeNames)
    indicesIRed <- match(uProbesIRed, uProbeNames)
    indicesII   <- match(uProbesII, uProbeNames)
    indicesX    <- match(uProbesX, uProbeNames)
    indicesY    <- match(uProbesY, uProbeNames)

    indList <- list(indicesIGrn, indicesIRed, indicesII, indicesX, indicesY)
    names(indList) <- c("IGrn", "IRed", "II","X","Y")


    ## Extraction of the quantiles
    mQuantiles               <- vector("list",5)
    betaQuantiles            <- vector("list", 5)
    methQuantiles            <- vector("list", 5)
    unmethQuantiles          <- vector("list", 5)
    cnQuantiles              <- vector("list", 5)
    names(mQuantiles)        <- c("IGrn", "IRed", "II","X","Y")
    names(betaQuantiles)     <- c("IGrn", "IRed", "II","X","Y")
    names(methQuantiles)     <- c("IGrn", "IRed", "II","X","Y")
    names(unmethQuantiles)   <- c("IGrn", "IRed", "II","X","Y")
    names(cnQuantiles)       <- c("IGrn", "IRed", "II","X","Y")

    nq <- 500
    probs <- seq(0,1,1/(nq-1))

    for (i in 1:5) {
        mQuantiles[[i]]      <- extractQuantiles(m[indList[[i]],], probs = probs)
        betaQuantiles[[i]]   <- extractQuantiles(beta[indList[[i]],], probs = probs)
        methQuantiles[[i]]   <- extractQuantiles(meth[indList[[i]],], probs = probs)
        unmethQuantiles[[i]] <- extractQuantiles(unmeth[indList[[i]],], probs = probs)
        cnQuantiles[[i]]     <- extractQuantiles(cn[indList[[i]],], probs = probs)
    }

    medianXU <-  unmethQuantiles$X[250,]
    medianXM <-  methQuantiles$X[250,]
    medianYU <-  unmethQuantiles$Y[250,]
    medianYM <-  methQuantiles$Y[250,]

    XYMedians <- list(medianXU = medianXU,
                      medianXM = medianXM,
                      medianYU = medianYU,
                      medianYM = medianYM
                      )

    return(list(
        mQuantiles = mQuantiles,
        betaQuantiles = betaQuantiles,
        methQuantiles = methQuantiles,
        unmethQuantiles = unmethQuantiles,
        cnQuantiles = cnQuantiles,
        greenControls = greenControls,
        redControls = redControls,
        XYMedians = XYMedians,
        oob = oob))
}
