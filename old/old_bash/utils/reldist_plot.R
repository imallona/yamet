#!/bin/env R


HOME <- '/home/imallona'
TASK="cg_shadows"
WD = file.path(HOME, 'mnt', 'nfs', TASK, 'data', 'ENCFF857QML')

setwd(WD)

## ctcfs <- c('ENCFF302QQX.bed', 'ENCFF559LUB.bed', 'ENCFF654IUM.bed', 'ENCFF719EWH.bed')


plot_reldist <- function(fn, tag_a, tag_b) {

 png( sprintf('%s_reldistplot.png', fn),
        width = 480, height = 480)
    par(mfrow = c(1,1),
        oma = c(5,4,1,1) + 0.1,
        mar = c(2,1,2,1) + 0.1,
        cex.main = 1.2,
        cex.lab = 1.2,
        cex.axis = 1.2)
    
    d <- read.table(fn,
                    sep = "\t", header = TRUE)

    plot(d$reldist, d$fraction, type = 'b',
         main = sprintf('%s vs\n %s', tag_a, tag_b),
         axes = TRUE,
         ylim = c(0, 1))
    

    title(xlab = sprintf("reldist %s vs\n %s", tag_a, tag_b),
          ylab = "frequency",
          outer = TRUE, line = 1)
    
    dev.off()
}


## ctcfs <- c('ENCFF302QQX.bed', 'ENCFF559LUB.bed', 'ENCFF654IUM.bed', 'ENCFF719EWH.bed')

## for (ctcf in ctcfs) {
##     for (item in c('plus', 'minus')) {
##         fn <- sprintf('%s_%s.reldist', ctcf, item)
##         plot_reldist(fn,  item, sprintf('ctcf %s', ctcf))
##     }
## }


## dnases <- c('ENCFF071XTK.bed', 'ENCFF421NEH.bed', 'ENCFF330CTT.bed', 'ENCFF006EIZ.bed',
##             'ENCFF722LRQ.bed', 'ENCFF257QKM.bed', 'ENCFF756XIR.bed', 'ENCFF113TVH.bed')

## for (dnase in dnases) {
##     for (item in c('plus', 'minus')) {
##         fn <- sprintf('%s_%s.reldist', dnase, item)
##         plot_reldist(fn,  item, sprintf('dnase %s', dnase))
##     }
## }

setwd(file.path(WD, 'reldist'))

items <- list.files(file.path(WD, 'reldist'), pattern = "*reldist$")

for (item in items) {
    tag <- gsub('.bed.reldist', '', item)
    tag_a <- strsplit(item, '__')[[1]][1]
    tag_b <- strsplit(item, '__')[[1]][3]
    
    plot_reldist(item, tag_a = tag_a, tag_b = tag_b)
}
