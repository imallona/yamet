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

    
    if(nrow(stretchesTable)>0){
      stretchesTable <- stretchesTable[,.(length_stretch=.N),
                                     by=c("temp_id", "bin")]}
    else{
      stretchesTable <- data.table(length_stretch=0)
    }
    
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
      lenMissing <- .getStretchLength(data, 1e4, naLength=TRUE)
      lenMissing <- lapply(lenMissing, function(l) pmax(l-1,0))

      lenCov <- .getStretchLength(data, 1e4, naLength=FALSE)
      nStr <- unlist(lapply(lenCov, length))
      lenCov <- lapply(lenCov, function(l) pmax(l-1,0))
      
      sampCells <- sample(which(nStr>20), nCells)
      lenCov <- lenCov[sampCells]
      lenMissing <- lenMissing[sampCells]
      
      paramMissing <- lapply(lenMissing, fitdistr, densfun="negative binomial")
      paramCov <- lapply(lenCov, fitdistr, densfun="negative binomial")
      
      paramMissing <- lapply(paramMissing, function(p){
        p <- p$estimate
        size <- p[1]
        prob <- size/(size+p[2])
        c(size, prob)})
      paramCov <- lapply(paramCov, function(p){
        p <- p$estimate
        size <- p[1]
        prob <- size/(size+p[2])
        c(size, prob)
      })
    }
    else{
      if(probParam<1){
        probCov <-  pmax(rep(probParam, nCells)+rnorm(nCells, sd=0.01), 0.01)
      }
      else{
        probCov <-  rep(1, nCells)
      }
      sizeCov <- rep(sizeParams[1], nCells)
      probMissing <- 1-probCov
      sizeMissing <- rep(sizeParams[2], nCells)
      paramCov <- lapply(1:length(probCov), function(i) c(sizeCov[i],
                                                          probMissing[i]))
     
      if(probParam<1){
        paramMissing <- lapply(1:length(probCov), function(i) c(sizeMissing[i],
                                                                probCov[i]))}
      else{
        paramMissing <- NULL
      }
    }
    
    sampNB <- function(params, nCpGs){
      nb <- rnbinom(n=nCpGs,
                    size=params[1],
                    prob=params[2])
      nb+1}
    
    if(!is.null(paramMissing)){
      lenMissingSamp <- lapply(paramMissing, sampNB, nCpGs=nCpGs)
      lenCovSamp <- lapply(paramCov, sampNB, nCpGs=nCpGs)
    }
    else{
      lenMissingSamp <- as.list(replicate(nCells, rep(0, nCpGs), simplify=FALSE))
      lenCovSamp <- as.list(replicate(nCells, nCpGs, simplify=FALSE))
    }
    
    # loop over and construct the per-cell data
    lenSamp <- lapply(1:nCells, function(i){
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
  else if(mode=="rand"){
      # estimate covered

      if(estimateParams  & !is.null(data)){
        cellIds <- setdiff(colnames(data), c("chr", "pos", "bin"))
        probs <- colSums(!is.na(data[, cellIds, with=FALSE]))/nrow(data)

        if(nCells>length(cellIds)) beRep <- TRUE else beRep <- FALSE
        probs <- sample(probs, nCells, replace=FALSE)
      }
      else{
        probs <- probParam
        probs <- rep(probs, nCells)
      }

      sampRand <- function(p, nCpGs){
        p <- c(1-p,p)
        ss <- sample(c(as.numeric(NA),1), size=nCpGs, prob=p, replace=TRUE)
        lenMissingSamp <- .getStretchLength(data.table(seq_bern=ss,
                                                       pos=1:length(ss),
                                                       chr=rep(1, length(ss))),
                                           nCpGs, naLength=TRUE)

        lenCovSamp <- .getStretchLength(data.table(seq_bern=ss,
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

  set.seed(seed)
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
  nCells <- length(lenCovSamp)
  cellSeeds <- sample.int(nCells*100, nCells)

  metPattern <- mapply(function(str, cellSeed){
    set.seed(cellSeed, kind = "L'Ecuyer-CMRG")
    seq <- lapply(str, markovchainSequence, markovChain)
    seq
  },lenCovSamp, cellSeeds, SIMPLIFY=FALSE)
  
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

      if(j<length(str) | (length(str)==1 & str[[j]]==0)){
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
#'@param probParam Probability parameter (`prob`) for [stats::rnbinom] (or [sample] when `mode="rand"`) to sample length of covered stretches.
#' If `probParam=1` a sequence with no-missing data will be simulated.
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
                          mode=c("nb", "rand"),
                          probParam=0.5,
                          sizeParams=c(1, 1),
                          estimateCovParams=FALSE,
                          estimateTransMat=FALSE,
                          transMat=NULL,
                          metTable=NULL,
                          seed=43)
{
  set.seed(seed)
  mode <- match.arg(mode, choices=c("nb", "rand"))
  
  # states of markovchain
  states=c('0', '1')
  
  lenSamp <- .simLengths(nCells, nCpGs, data=metTable, mode=mode,
                         estimateParams=estimateCovParams,
                         probParam=probParam,
                         sizeParams=sizeParams,
                         seed=seed)

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
