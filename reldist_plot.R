#!/bin/env R


HOME <- '/home/imallona'
TASK="cg_shadows"
WD = file.path(HOME, TASK, 'data')

setwd(WD)

ctcfs <- c('ENCFF302QQX.bed', 'ENCFF559LUB.bed', 'ENCFF654IUM.bed', 'ENCFF719EWH.bed')

d <- list()

for (ctcf in ctcfs) {

    png(file.path(WD, sprintf('%s_plus_reldist.png', ctcf)),
        width = 480, height = 480)
    par(mfrow = c(1,1),
        oma = c(5,4,1,1) + 0.1,
        mar = c(2,1,2,1) + 0.1
        )
    
    d[[as.character(ctcf)]]<- read.table(sprintf('%s_plus.reldist', ctcf),
                                         sep = "\t", header = TRUE)


    ## ymax = max(unlist(sapply(d[[as.character(ctcf)]], function(x) return(x$fraction))))*1.1
    
    plot(d[[ctcf]]$reldist, d[[ctcf]]$fraction, type = 'b',
         main = sprintf('%s', ctcf),
         axes = TRUE,
         ylim = c(0, 0.2))
    

    title(xlab = sprintf("relative distance between plus agreements and ctcf %s", ctcf),
          ylab = "frequency",
          outer = TRUE, line = 1)
    
    dev.off()

}


for (ctcf in ctcfs) {

    png(file.path(WD, sprintf('%s_minus_reldist.png', ctcf)),
        width = 480, height = 480)
    par(mfrow = c(1,1),
        oma = c(5,4,1,1) + 0.1,
        mar = c(2,1,2,1) + 0.1
        )
    
    d[[as.character(ctcf)]]<- read.table(sprintf('%s_minus.reldist', ctcf),
                                         sep = "\t", header = TRUE)


    ## ymax = max(unlist(sapply(d[[as.character(ctcf)]], function(x) return(x$fraction))))*1.1
    
    plot(d[[ctcf]]$reldist, d[[ctcf]]$fraction, type = 'b',
         main = sprintf('%s', ctcf),
         axes = TRUE,
         ylim = c(0, 0.2))
    

    title(xlab = sprintf("relative distance between minus agreements and ctcf %s", ctcf),
          ylab = "frequency",
          outer = TRUE, line = 1)
    
    dev.off()

}

dnases <- c('ENCFF071XTK.bed', 'ENCFF421NEH.bed', 'ENCFF330CTT.bed', 'ENCFF006EIZ.bed',
            'ENCFF722LRQ.bed', 'ENCFF257QKM.bed', 'ENCFF756XIR.bed', 'ENCFF113TVH.bed')


d <- list()

for (dnase in dnases) {

    png(file.path(WD, sprintf('%s_plus_reldist.png', dnase)),
        width = 480, height = 480)
    par(mfrow = c(1,1),
        oma = c(5,4,1,1) + 0.1,
        mar = c(2,1,2,1) + 0.1
        )
    
    d[[as.character(dnase)]]<- read.table(sprintf('%s_plus.reldist', dnase),
                                         sep = "\t", header = TRUE)


    ## ymax = max(unlist(sapply(d[[as.character(dnase)]], function(x) return(x$fraction))))*1.1
    
    plot(d[[dnase]]$reldist, d[[dnase]]$fraction, type = 'b',
         main = sprintf('%s', dnase),
         axes = TRUE,
         ylim = c(0, 0.2))
    

    title(xlab = sprintf("relative distance between plus agreements and dnase %s", dnase),
          ylab = "frequency",
          outer = TRUE, line = 1)
    
    dev.off()

}


for (dnase in dnases) {

    png(file.path(WD, sprintf('%s_minus_reldist.png', dnase)),
        width = 480, height = 480)
    par(mfrow = c(1,1),
        oma = c(5,4,1,1) + 0.1,
        mar = c(2,1,2,1) + 0.1
        )
    
    d[[as.character(dnase)]]<- read.table(sprintf('%s_minus.reldist', dnase),
                                         sep = "\t", header = TRUE)


    ## ymax = max(unlist(sapply(d[[as.character(dnase)]], function(x) return(x$fraction))))*1.1
    
    plot(d[[dnase]]$reldist, d[[dnase]]$fraction, type = 'b',
         main = sprintf('%s', dnase),
         axes = TRUE,
         ylim = c(0, 0.2))
    

    title(xlab = sprintf("relative distance between minus agreements and dnase %s", dnase),
          ylab = "frequency",
          outer = TRUE, line = 1)
    
    dev.off()

}
