# Functionalities for simulating methylation patterns with differing coverage
#'@author: Emanuel Sonder

.getStretchLength <- function(metTable,
                              nCpGs=NULL,
                              subsetRange=NULL,
                              naLength=TRUE){
  cellIds <- setdiff(names(metTable), c("pos", "chr", "bin"))

  if(!is.null(subsetRange))
  {
    metTable <- metTable[subsetRange[1]:subsetRange[2],]
  }

  if(!is.null(nCpGs))
  {
    metTable[, bin:=cut(pos, seq(min(pos), max(pos)+nCpGs, nCpGs),
                        include.lowest=TRUE), by=chr]
  }
  else{
    metTable$bin <- 1
  }

  setorder(metTable, chr, pos)
  pos <- metTable$pos
  bin <- metTable$bin

  callLength <- function(col, pos, naLength=TRUE){

    # mark stretches
    if(naLength){
      isEnd <- fifelse(!is.na(data.table::shift(col, n=1, type="lead")) &
                         is.na(col), TRUE, FALSE)
    }
    else{
      isEnd <- fifelse(is.na(data.table::shift(col, n=1, type="lead")) &
                         !is.na(col), TRUE, FALSE)
    }


    tempId <- fifelse(isEnd, pos, as.integer(NA))
    tempId <- nafill(tempId, type="nocb")

    stretchesTable <- data.table(temp_id=tempId,
                                 bin=bin,
                                 rate=col)

    # get length of NA stretches
    if(naLength){
      stretchesTable <- subset(stretchesTable, is.na(rate))
    }
    else{
      stretchesTable <- subset(stretchesTable, !is.na(rate))
    }


    stretchesTable <- stretchesTable[,.(length_stretch=.N),
                                     by=c("temp_id", "bin")]

    return(stretchesTable$length_stretch)
  }

  stretchesLength <- lapply(metTable[,cellIds, with=FALSE], callLength, pos=metTable$pos, naLength)
  #stretchesLength <- unlist(stretchesLength)

  return(stretchesLength)
}

.simLengths <- function(nCells, nCpGs, data=NULL, mode=c("nb", "random"),
                        estimateParams=FALSE,
                        probParam=NULL,
                        sizeParams=NULL,
                        seed=43){
  set.seed(seed)
  if(mode=="nb"){
    if(estimateParams & !is.null(data)){
      lenMissing <- .getStretchLength(data, nCpGs, naLength=TRUE)
      sampCells <- sample(1:length(lenMissing), nCells)
      lenMissing <- lenMissing[sampCells]
      lenMissing <- lapply(lenMissing, function(l) l-1)

      lenCov <- .getStretchLength(data, nCpGs, naLength=FALSE)
      lenCov <- lenCov[sampCells]
      lenCov <- lapply(lenCov, function(l) l-1)

      paramMissing <- lapply(lenMissing, fitdistr, densfun="negative binomial")
      paramCov <- lapply(lenCov, fitdistr, densfun="negative binomial")

      sampNB <- function(params, nCpGs){
        nb <- rnbinom(n=nCpGs,
                      size=params$estimate[1],
                      mu=params$estimate[2])
        nb+1}
    }
    else{
      probCov <-  pmin(rep(probParam, nCells)+rnorm(100, sd=0.01), 0.99)
      sizeCov <- pmax(rep(sizeParams[1], nCells)+rnorm(100, sd=0.5),0.5)
      paramCov <- lapply(1:length(probCov), function(i) c(sizeCov[i],
                                                          probCov[i]))
      
      probMissing <- 1-probCov
      sizeMissing <- pmax(rep(sizeParams[2], nCells)+rnorm(100, sd=0.5),0.5)
      paramMissing <- lapply(1:length(probCov), function(i) c(sizeMissing[i],
                                                              probMissing[i]))
  
      sampNB<- function(params, nCpGs){
        nb <- rnbinom(n=nCpGs,
                      size=params[1],
                      prob=params[2])
        nb+1}
    }
    lenMissingSamp <- lapply(paramMissing, sampNB, nCpGs=nCpGs)
    lenCovSamp <- lapply(paramCov, sampNB, nCpGs=nCpGs)

    # loop over and construct the per-cell data
    lenSamp <- lapply(1:nCells, function(i){
      lenMissingSamp
      strDt <- data.table(len_miss=lenMissingSamp[[i]],
                          len_cov=lenCovSamp[[i]])

      strDt[,sum_miss:=cumsum(len_miss)]
      strDt[,sum_cov:=cumsum(len_cov)]
      strDt[,sum_tot:=sum_miss+sum_cov]
      strDt <- subset(strDt, sum_tot<=nCpGs)
      nDiff <- nCpGs-max(strDt$sum_tot)

    return(list("len_missing"=c(strDt$len_miss, nDiff),
                "len_cov"=strDt$len_cov))
    })
  }
  else if(mode=="random"){
      # estimate covered

      if(estimateParams  & !is.null(data)){
        cellIds <- setdiff(colnames(metTable), c("chr", "pos", "bin"))
        probs <- colSums(!is.na(metTable[, cellIds, with=FALSE]))/nrow(metTable)

        if(nCells>length(cellIds)) beRep <- TRUE else beRep <- FALSE
        probs <- sample(probs, nCells, replace=FALSE)
      }
      else{
        probs <- probParam[1]
        probs <- rep(probs, nCells)
      }

      sampRand <- function(p, nCpGs){
        p <- c(1-p,p)
        ss <- sample(c(as.numeric(NA),1), size=nCpGs, prob=p, replace=TRUE)
        lenMissingSamp <- .getStretchLength(data.table(a=ss,
                                                       pos=1:length(ss),
                                                       chr=rep(1, length(ss))),
                                           nCpGs, naLength=TRUE)
        lenCovSamp <- .getStretchLength(data.table(a=ss,
                                                   pos=1:length(ss),
                                                   chr=rep(1, length(ss))),
                                       nCpGs, naLength=FALSE)
        return(list(len_missing=unlist(lenMissingSamp),
                    len_cov=unlist(lenCovSamp)))
        }
      lenSamp <- lapply(probs, sampRand, nCpGs=nCpGs)
  }

  lenCovSamp <- lapply(lenSamp, function(samp) samp$len_cov)
  lenMissingSamp <- lapply(lenSamp, function(samp) samp$len_missing)

  return(list("length_covered"=lenCovSamp,
              "length_missing"=lenMissingSamp))
}

.simCovData <- function(lenCovSamp,
                        transMat=matrix(rep(0.5, 4), nrow=2, ncol=2),
                        estimateTransMat=FALSE,
                        data=NULL,
                        states=c('0', '1'),
                        seed=43){

  if(estimateTransMat & !is.null(data)){
    cellIds <- setdiff(colnames(data), c("chr", "pos", "bin"))
    seqData <- lapply(cellIds, function(cell){data[[cell]]})
    seqData <- lapply(seqData, function(cell){as.character(cell)})
    #seqData <- lapply(seqData, function(cell){cell[!(cell %in% states)] <- NA
    #cell})
    empEst <- markovchainFit(seqData, method="mle", possibleStates=states)
    transMat <- empEst$estimate@transitionMatrix
  }

  markovChain <- new('markovchain',
                     transitionMatrix=transMat,
                     states=states)

  # generate markov sequence
  set.seed(seed, kind = "L'Ecuyer-CMRG")

  metPattern <- lapply(lenCovSamp, function(str){
    seq <- lapply(str, markovchainSequence, markovChain)
    seq
  })

  metPattern <- lapply(metPattern, function(str){
    seq <- lapply(str, as.numeric)
    seq})

  return(metPattern)
}

.simData <- function(lenCovSamp, lenMissingSamp, nCells,
                     transMat=NULL,
                     estimateTransMat=FALSE,
                     data=NULL,
                     states=c('0', '1'),
                     seed=43){

  set.seed(seed)
  # simulate covered data
  metPattern <- .simCovData(lenCovSamp,
                            transMat=transMat,
                            estimateTransMat=estimateTransMat,
                            data=data,
                            states=states,
                            seed=seed)

  # put together with missing stretches
  simCells <- lapply(1:nCells, function(i){

    str <- lenMissingSamp[[i]]
    seq <- lapply(1:length(str), function(j){

      if(j<length(str)){
        covSeq <- metPattern[[i]][[j]]
        }
      else{
        covSeq <- NULL
      }
      c(rep(NA, str[[j]]), covSeq)
    })
    simCell <- data.table(unlist(seq))
    colnames(simCell) <- paste("sim_cell", i, sep="_")
    simCell
  })

  simDat <-  Reduce(cbind,simCells[-1], simCells[[1]])
}

#'@title simMetPattern
#'@description one feature at a time
#'@param nCpGs length of regions in number of CpGs to simulate
#'@param nCells number of cells to simulate
#'@param mode Mode for simulating lengths of missing and covered stretches,
#'either "nb" (sampling from a negative binomial) or "rand" (repeated bernoulli)
#'@param probParam Probability parameter (`prob`) for [stats::rnbinom] to sample length of covered stretches.
#'@param sizeParams vector of size parameters for rnbinom to sample length of stretches.
#' First element corresponds to the size parameter for sampling covered stretches, second for missing stretches (with `prob=1-probParam`).
#'@param states states to draw from
#'@param estimateCovParams Should coverage parameters be estimated from provided methylation data
#'@param estimateTransMat Should transition matrix be estimated from provided methylation data
#'@param transMat 2 by 2 transition matrix for markov chain
#'@param metTable methylation data to estimate parameters from
#'@param seed random seed
#'@return
simMetPattern <- function(nCpGs, nCells,
                          mode=c("nb", "random"),
                          probParam=0.5,
                          sizeParams=c(1, 1),
                          estimateCovParams=FALSE,
                          estimateTransMat=FALSE,
                          transMat=NULL,
                          metTable=NULL,
                          seed=43)
{
  set.seed(seed)
  mode <- match.arg(mode, choices=c("nb", "random"))
  
  # states of markovchain
  states=c('0', '1')
  
  lenSamp <- .simLengths(nCells, nCpGs, data=metTable, mode=mode,
                         estimateParams=estimateCovParams,
                         probParam=probParam,
                         sizeParams=sizeParams,
                         seed=seed)
  # quick fix
  lenMissing <- lapply(lenSamp$length_missing, function(str){
    str <- fifelse(!is.finite(str),nCpGs, str)
  })              

  metTable <- .simData(lenSamp$length_covered,
                       lenMissing,
                       nCells=nCells,
                       transMat=transMat,
                       estimateTransMat=estimateTransMat,
                       data=metTable,
                       states=states,
                       seed=seed)

  metTable$pos <- seq(1,nCpGs,1)
  metTable$chr <- "chrSim"

  return(metTable)
}
