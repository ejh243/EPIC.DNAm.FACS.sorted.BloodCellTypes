#' Function to load reference data included in this package and process with user supplied DNA methylation data to estimate cell composition for blood. This is an implementaion of the Houseman et al (2012) regression calibration approach algorithm for deconvoluting heterogeneous tissue sources like blood. This function is adapted from the minfi package. This function will take an RGChannelSet from a DNA methylation (DNAm) study of blood, and return the relative proportions of CD4+ and CD8+ T-cells, monocytes, granulocytes, and b-cells in each sample. The function requires an RGChannelSet as prior to the estimation step, it normalises the referecne data and user supplied data together. Although the reference data was generated with the Illumina EPIC microarray, users can provde 450K datasets for estimation. 

#' 

#' @param userData is an RGChannelSet
#' @param EPIC is logical indicate whether userData were profiled with the 450K array
#' @param processMethod Specify how the user and reference data should be processed together. Default input "auto" will use preprocessQuantile in line with the existing literature. Set it to the name of a preprocessing function as a character if you want to override it, like "preprocessFunnorm". See minfi for available options.
#' @param probeSelect Specify how probes should be selected to distinguish cell types. Options include "both", which selects an equal number (50) of probes (with F-stat p-value < 1E-8) with the greatest magnitude of effect from the hyper- and hypo-methylated sides, and "any", which selects the 100 probes (with F-stat p-value < 1E-8) with the greatest magnitude of difference regardless of direction of effect. Default input "auto" will use "both" for blood, in line with previous versions of this function and/or our recommendations.
#' @param cellTypes A vector to specify which cell types the user wants to estimate. A subset of c("CD8T","CD4T", "Bcell","Mono","Gran") is allowed.
#' @param returnAll Logical inicating whether the composition table and the normalized user supplied data be returned.
#' @param meanPlot Logical indicating whether to plot the average DNA methylation across the cell-type discrimating probes within the mixed and sorted samples. This should be used to check for large batch effects in the data, reducing the confidence placed in the composition estimates. 
#' @param verbose Should the function be verbose?
#' @return A matrix containing the estimated proportion of each cell type for each sample. Columns contain cell types while rows contain samples. 
#' @return If returnAll=TRUE a list of a count matrix (see previous paragraph), a composition table and the normalized user data in form of a GenomicMethylSet.
#' @return If meanPlot=TRUE A plot depicting the average DNA methylation across the cell-type discrimating probes in both the provided and sorted data is produced. The means from the provided heterogeneous samples should be within the range of the sorted samples. If the sample means fall outside the range of the sorted means, the cell type estimates will inflated to the closest cell type. Note that we quantile normalize the sorted data with the provided data to reduce these batch effects.
#' @examples
#' \dontrun{ 
#' counts<-estimateCellCountsEPIC(RGSet, EPIC=TRUE, cellTypes=c("Bcells", "CD4Tcells", "CD8Tcells", "Monocytes", "Granulocytes"))
#' }


estimateCellCountsEPIC<-function(userData, EPIC = TRUE, processMethod = "auto",probeSelect = "auto", cellTypes = c("Bcells","CD4Tcells","CD8Tcells","Monocytes","Granulocytes"), returnAll = FALSE, meanPlot = FALSE, verbose = TRUE){

	if ((processMethod == "auto"))
        processMethod <- "preprocessQuantile"
    processMethod <- get(processMethod)
     if ((probeSelect == "auto")){
        probeSelect <- "both"}
	if(EPIC){
		referenceRGset<-RGSet.ref
	} else {
		referenceRGset<-convertArray(RGSet.ref,outType = "IlluminaHumanMethylation450k",verbose = FALSE)
	}
	if(verbose) message("Combining user data with reference (flow sorted) data.\n")
    newpd <- DataFrame(sampleNames = c(colnames(userData), colnames(referenceRGset)),
                       studyIndex = rep(c("user", "reference"),
                                        times = c(ncol(userData), ncol(referenceRGset))))
    referencePd <- pData(referenceRGset)
	if(EPIC){
		combinedRGset <- combineArrays(userData, referenceRGset, outType = "IlluminaHumanMethylationEPIC")
	
	} else {
		combinedRGset <- combineArrays(userData, referenceRGset, outType = "IlluminaHumanMethylation450k")
	}
    pData(combinedRGset) <- newpd
    colnames(combinedRGset) <- newpd$sampleNames
    rm(referenceRGset)
    
    if(verbose) message("Processing user and reference data together.\n")
    
    combinedMset <- processMethod(combinedRGset) 
    rm(combinedRGset)
    
    ## Extracts normalized reference data 
    referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
    colData(referenceMset) <- as(referencePd, "DataFrame")
    mSet <- combinedMset[, combinedMset$studyIndex == "user"]
    colData(mSet) <- as(pData(userData), "DataFrame")
    rm(combinedMset)

	if(verbose) cat("Picking probes for composition estimation.\n")
	
	compData <- pickCompProbes(referenceMset, cellTypes = cellTypes)
	coefs <- compData$coefEsts
		
	if(verbose) cat("Estimating composition.\n")
	counts <- projectCellType(getBeta(userData)[rownames(coefs), ], coefs)
	rownames(counts) <- sampleNames(userData)
	
	if (meanPlot) {
        smeans <- compData$sampleMeans
        smeans <- smeans[order(names(smeans))]
        sampleMeans <- c(colMeans(getBeta(mSet)[rownames(coefs), ]), smeans)
        
        sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
        plot(sampleMeans, pch = 21, bg = sampleColors)
        legend("bottomleft", c("blood", levels(factor(names(smeans)))),
               col = 1:7, pch = 15)
    }
    if(returnAll) {
        list(counts = counts, compTable = compData$compTable,
             normalizedData = mSet)
    } else {
        counts
    }
}
validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]
    
    if(is.null(L.forFstat)) {
        L.forFstat <- diag(sizeModel)[-1,] # All non-intercept coefficients
        colnames(L.forFstat) <- colnames(xTest) 
        rownames(L.forFstat) <- colnames(xTest)[-1] 
    }

    ## Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()

    if(verbose)
        cat("[validationCellType] ")
    for(j in 1:M) { # For each CpG
        ## Remove missing methylation values
        ii <- !is.na(Y[j,])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j,]
        
        if(j%%round(M/10)==0 && verbose)
            cat(".") # Report progress
        
        try({ # Try to fit a mixed model to adjust for plate
            if(!is.null(modelBatch)) {
                fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
                OLS <- inherits(fit,"try-error") # If LME can't be fit, just use OLS
            } else
                OLS <- TRUE

            if(OLS) {
                fit <- lm(modelFix, data=pheno[ii,])
                fitCoef <- fit$coef
                sigmaResid[j] <- summary(fit)$sigma
                sigmaIcept[j] <- 0
                nClusters[j] <- 0
            } else { 
                fitCoef <- fit$coef$fixed
                sigmaResid[j] <- fit$sigma
                sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
                nClusters[j] <- length(fit$coef$random[[1]])
            }
            coefEsts[j,] <- fitCoef
            coefVcovs[[j]] <- vcov(fit)
            
            useCoef <- L.forFstat %*% fitCoef
            useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
            Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
        })
    }
    if(verbose)
        cat(" done\n")
    ## Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) <- rownames(Y)
    colnames(coefEsts) <- names(fitCoef)
    degFree <- nObserved - nClusters - sizeModel + 1

    ## Get P values corresponding to F statistics
    Pval <- 1-pf(Fstat, sizeModel, degFree)
    
    out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
                sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval,
                orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, nObserved=nObserved,
                degFree=degFree)
    
    out
}

pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50, compositeCellType = "Blood", probeSelect = "auto") {
    splitit <- function(x) {
        split(seq(along=x), x)
    }
    
    p <- getBeta(mSet)
    pd <- as.data.frame(colData(mSet))
    if(!is.null(cellTypes)) {
        if(!all(cellTypes %in% pd$CellType))
            stop("elements of argument 'cellTypes' is not part of 'mSet$CellType'")
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep,]
        p <- p[,keep]
    }
    ## make cell type a factor 
    pd$CellType <- factor(pd$CellType)
    ffComp <- rowFtests(p, pd$CellType)
    prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
    r <- matrixStats::rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range") 
    tIndexes <- splitit(pd$CellType)
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })
    
    if (probeSelect == "any"){
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]      
            c(rownames(yAny)[1:(numProbes*2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yUp <- y[order(y[,"dm"], decreasing=TRUE),]
            yDown <- y[order(y[,"dm"], decreasing=FALSE),]
            c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
        })
    }
    
    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes,]
    
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse="+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType-1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if(ncol(phenoDF) == 2) { # two group solution
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
    } else { # > 2 group solution
        tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
        coefEsts <- tmp$coefEsts
    }
    
    out <- list(coefEsts = coefEsts, compTable = compTable,
                sampleMeans = pMeans)
    return(out)
}

projectCellType <- function(Y, coefCellType, contrastCellType=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
    if(is.null(contrastCellType))
        Xmat <- coefCellType
    else
        Xmat <- tcrossprod(coefCellType, contrastCellType) 
    
    nCol <- dim(Xmat)[2]
    if(nCol == 2) {
        Dmat <- crossprod(Xmat)
        mixCoef <- t(apply(Y, 2, function(x) { solve(Dmat, crossprod(Xmat, x)) }))
        colnames(mixCoef) <- colnames(Xmat)
        return(mixCoef)
    } else {
        nSubj <- dim(Y)[2]
        
        mixCoef <- matrix(0, nSubj, nCol)
        rownames(mixCoef) <- colnames(Y)
        colnames(mixCoef) <- colnames(Xmat)
        
        if(nonnegative){
            if(lessThanOne) {
                Amat <- cbind(rep(-1, nCol), diag(nCol))
                b0vec <- c(-1, rep(0, nCol))
            } else {
                Amat <- diag(nCol)
                b0vec <- rep(0, nCol)
            }
            for(i in 1:nSubj) {
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), Amat, b0vec)$sol
            }
        } else {
            for(i in 1:nSubj) {
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i])
            }
        }
        return(mixCoef)
    }
}

library(minfi)
library(genefilter)
library(quadprog)



