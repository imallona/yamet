#!/usr/bin/env R
##
## Cometh agreement scoring simulations
##
## Izaskun Mallona
## 12th nov 2018

library('latex2exp')
library(entropy)
library(Cairo)

NROWS <- 1e6
WD <- file.path('/home/imallona', 'cg_shadows', 'data', 'ENCFF857QML')
meth_fn <- file.path(WD, 'mini.CG.2_meth_colored.bed')
entropy_fn <- file.path(WD, 'mini.CG.2_entropy_colored.bed')

beta2m <- function(beta) {
    m <- log2(beta/(1 - beta))
    m
}


d <- list (entropy = read.table(entropy_fn,
                                colClasses = c(rep("NULL", 3), "factor", "numeric",
                                               "NULL", "factor"),
                                nrows = NROWS),
           meth = read.table(meth_fn,
                             colClasses = c(rep("NULL", 4), "numeric", "NULL", "factor"),
                             nrows = NROWS))

colors <- read.table(text = "Art
Ctcf
CtcfO
DnaseD
DnaseU
Elon
ElonW
Enh
EnhF
EnhW
EnhWF
FaireW
Gen3'
Gen5'
H4K20
Low
Pol2
PromF
PromP
Quies
Repr
ReprD
ReprW
Tss
TssF")$V1


                                
colnames(d$entropy) <- c('code', 'entropy', 'hmm')
colnames(d$meth) <- c('beta', 'hmm')
d$meth$beta <- d$meth$beta/1000

## d$entropy$mm <- sapply(strsplit(as.character(d$entropy$code), ';'),
##                        function(x) return(gsub('MM', '', x)))


tmp  <- sapply(strsplit(as.character(d$entropy$code), ';'),
               function(x) return(substr(x,3, nchar(x))))

d$entropy$mm <- as.numeric(tmp[1,])
d$entropy$mu <- as.numeric(tmp[2,])
d$entropy$um <- as.numeric(tmp[3,])
d$entropy$uu <- as.numeric(tmp[4,])

rm(tmp)

for (color in colors) {
    d$entropy[,color] <- NA
    d$entropy[,color] <- grepl(color, as.character(d$entropy$hmm))
}




d$entropy$in_boundary <- apply(d$entropy[,as.character(colors)], 1, function(x) sum(x) > 1)
table(d$entropy$in_boundary)

## wide to long?

## summary(d$entropy[!d$entropy$in_boundary, 'entropy'])

## fit1 <- aov(entropy~hmm, data = d$entropy)
## summary(fit1)

stripchart(d$entropy$entropy~ d$entropy$hmm, jitter = 0.1,
           vertical = TRUE,
           pch = 21)


myFun <- function(x) {
  c(min = min(x), max = max(x), 
    mean = mean(x), median = median(x), 
    std = sd(x))
}

tapply(d$entropy$entropy, d$entropy$hmm, function(x) quantile(x, probs = 0.99))

tapply(d$entropy$entropy, d$entropy$hmm, function(x) var(x))

##########

summary(d)
d$coverage <- d$MM + d$UM + d$MU + d$UU
d$beta <- (d$MM + 0.5*d$UM + 0.5*d$MU) / d$coverage
d$m <- beta2m(d$beta)

d$score1 <- (d$MM + d$UU) / d$coverage
d$score2 <- (d$MM) / d$coverage

d$discordant <- d$MU + d$UM
## d$score3 <- (2*d$MM + d$UU) / d$coverage
d$entropy_four <- apply(d[,c('MM','MU', 'UM', 'UU')], 1, entropy)
d$entropy_three <- apply(d[,c('MM','discordant', 'UU')], 1, entropy)

d$distance <- d$pos2 - d$pos1
d$log10dist <- log10(d$distance)
head(d$entropy)

## targets <- c('coverage', 'distance', 'beta',
##              'm', 'score1', 'score2',
##              'entropy_four', 'entropy_three')

targets <- c('coverage',
             'distance', 'log10dist',
             'beta',
             'm',
             'score1',
             ## 'score2',
             'entropy_four')#,
             ## 'entropy_three')

png(file.path(WD, 'exploratory_1.png'), width = 700,
         height = 700)

plot(d[,targets], pch = 20)#, col = rgb(0, 0, 1, 0.5))
dev.off()





    
pdf(file.path(WD, 'exploratory_2.pdf'), width = 14)


## scatterplotMatrix(d[,targets], main="", diagonal = 'density', plot.points = FALSE)

NPLOT = 150


## score3


## plot(y = jitter(d$score3[1:NPLOT]),
##      x = jitter(d$beta[1:NPLOT]),
##      main = TeX('$score3= \\frac{2*MM + UU}{cov} vs \\beta$'),
##      xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
##      ylab = 'score3 with jitter') 


plot(y = d$score3[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 2,
     ylab= 'score',
     xlab = 'tuple index (sorted by coordinates)',
     xlim = c(0, NPLOT + 0.2*NPLOT))

lines(y = d$beta[1:NPLOT], x= 1:NPLOT, type = 'l', col = 'blue', lwd = 2)
lines(y = log10(d$coverage[1:NPLOT]), x= 1:NPLOT, type = 'p', lwd = 1, col = 'darkred',
      pch = 20)

legend("topright",
       legend = c(TeX('$score3= \\frac{2*MM + UU}{cov}$'),
                  TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$'),
                  'log10(number of reads)'),
       col = c('black', 'blue', 'darkred'),
       lwd = c(2,2,1),
       lty = c(1,1,0),
       pch = c(NA, NA, 20), 
       ## xjust = 1, yjust = 1,
       title = "Scores")


## jittering
plot(y = jitter(d$score3[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     ylab = 'score3 with jitter',
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     main =  TeX('$score3= \\frac{2*MM + UU}{cov}$'))




## on score1
plot(y = jitter(d$score1[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     main = TeX('$score1= \\frac{MM + UU}{cov} vs \\beta$'),
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     ylab = 'score1 with jitter') 


plot(y = d$score1[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 2,
     ylab= 'score1',
     xlab = 'tuple index (sorted by coordinates)',
     xlim = c(0, NPLOT + 0.2*NPLOT))

lines(y = d$beta[1:NPLOT], x= 1:NPLOT, type = 'l', col = 'blue', lwd = 2)
lines(y = log10(d$coverage[1:NPLOT]), x= 1:NPLOT, type = 'p', lwd = 1, col = 'darkred',
      pch = 20)

legend("topright",
       legend = c(TeX('$score1= \\frac{MM + UU}{cov}$'),
                  TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$'),
                  'log10(number of reads)'),
       col = c('black', 'blue', 'darkred'),
       lwd = c(2,2,1),
       lty = c(1,1,0),
       pch = c(NA, NA, 20), 
       ## xjust = 1, yjust = 1,
       title = "Scores")


## jittering
plot(y = jitter(d$score1[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     ylab = 'score1 with jitter',
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     main = TeX('$score1= \\frac{MM + UU}{cov}$'))





## on score2
plot(y = jitter(d$score2[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     main = TeX('$score2= \\frac{MM}{cov} vs \\beta$'),
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     ylab = 'score2 with jitter') 


plot(y = d$score2[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 2,
     ylab= 'score2',
     xlab = 'tuple index (sorted by coordinates)',
     xlim = c(0, NPLOT + 0.2*NPLOT))

lines(y = d$beta[1:NPLOT], x= 1:NPLOT, type = 'l', col = 'blue', lwd = 2)
lines(y = log10(d$coverage[1:NPLOT]), x= 1:NPLOT, type = 'p', lwd = 1, col = 'darkred',
      pch = 20)

legend("topright",
       legend = c(TeX('$score2= \\frac{MM}{cov}$'),
                  TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$'),
                  'log10(number of reads)'),
       col = c('black', 'blue', 'darkred'),
       lwd = c(2,2,1),
       lty = c(1,1,0),
       pch = c(NA, NA, 20), 
       ## xjust = 1, yjust = 1,
       title = "Scores")


## jittering
plot(y = jitter(d$score2[1:NPLOT]),
     x = jitter(d$beta[1:NPLOT]),
     ylab = 'score2 with jitter',
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     main = TeX('$score2= \\frac{MM}{cov}$'))


dev.off()






## lineplots depicting score 1, 3 and beta



plot(y = d$score1[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 2,
     ylab= 'score1',
     xlab = 'tuple index (sorted by coordinates)',
     xlim = c(0, NPLOT + 0.2*NPLOT),
     ylim = c(0,1))

lines(y = d$beta[1:NPLOT], x= 1:NPLOT, type = 'l', col = 'blue', lwd = 2)
lines(y = d$score2[1:NPLOT], x= 1:NPLOT, type = 'l', lwd = 1, col = 'darkgreen')

legend("topright",
       legend = c(TeX('$score1= \\frac{MM + UU}{cov}$'),
                  TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$'),
                  '$score2 = \\frac{MM}{cov}$'),
       col = c('black', 'blue', 'darkgreen'),
       lwd = c(2,2,2),
       lty = c(1,1,1),
       ## xjust = 1, yjust = 1,
       title = "Scores")

     
## plot(y = jitter(d$score1[1:NPLOT] - d$beta[1:NPLOT]), x = jitter(d$beta[1:NPLOT]))


## plot(y = d$score1[1:NPLOT], x = 1:NPLOT, type = 'l')
## lines(y = d$beta[1:NPLOT], x = 1:NPLOT, type = 'l',  col = 'blue')
    
## highlighting dnameth


## jittering (not qvalue adjusted)
## checking the soundness of a dmr approach

png('exploratory_3.png')
d$dmr <- d$score1 - d$beta
d$dmr_color <- 'black'
d$dmr_color[d$dmr >= 0.1] <- 'blue'
d$dmr_color[d$dmr < -0.1] <- 'red'

plot(y = jitter(d$score1),
     x = jitter(d$beta),
     ylab = 'score1 with jitter',
     xlab = TeX('$ \\beta = \\frac{MM + .5*MU + .5*UM}{cov}$ with jitter'),
     main = TeX('$score1= \\frac{MM + UU}{cov}$'),
     col = d$dmr_color)


legend("bottomleft",
       legend = c('nochange', TeX('$score1 \\geq \\beta$'), TeX('$score1 < \\beta$')),
       col = c('black', 'blue', 'red'),
       pch = 1,
       ## xjust = 1, yjust = 1,
       title = "Scores")
dev.off()


## entropy measures

  
d$entropy <- apply(d[,c('MM','MU', 'UM', 'UU')], 2, entropy)
head(d$entropy)

plot(d[,c(targets, 'entropy')])

## set.seed(1)
## n = 500
## x = round(rbeta(n, 0.5, 0.2), 2)

## hist(x, prob = TRUE, main = 'beta')
## rug(x)
## curve(dbeta(x, .5, .4), lwd=2, col="blue", add=TRUE)

## ## cov <- floor(abs(rnorm(n,20,10)))
    
## ## fd <- data.frame(beta = x,
## ##                  cov = cov)
## ## fd$cov_meth <- floor(fd$beta * fd$cov)

## ## ## a third of the cases are MM tuples
## ## fd$mm <- fd$um <- fd$mu <- floor(fd$cov_meth/3)
## ## fd$uu <- fd$cov - (fd$mm + fd$um + fd$mu)

## ## summary(fd)

## ## hist(fd$cov_meth)
## ## plot(fd)



    
## ## fd <- data.frame(beta = x,
## ##                  mm = floor(abs(rnorm(n,5,2))),
## ##                  um = floor(abs(rnorm(n,5,2))),
## ##                  mu =  floor(abs(rnorm(n,5,2))))

## ## fd$meth <- (2*fd$mm) + fd$um + fd$mu
## ## fd$cov <- fd$meth / fd$beta
## ## fd$uu <- fd$cov - (fd$mm + fd$mu + fd$mm)

## ## fd <- fd[is.finite(fd$uu) & fd$uu <100,]

## ## plot(fd)

## ## ## the beta vs coverage is shitty


    
## fd <- data.frame(beta = x,
##                  cov = floor(abs(rnorm(n,25,3))),
##                  mm = floor(abs(rnorm(n,5,2))),
##                  uu =  floor(abs(rnorm(n,5,2))))

                 

## fd$missing_meth <- fd$cov* fd$beta
## fd$mu <- round(fd$missing_meth/2)
## fd$um <- round(fd$missing_meth/2)
    
## ## fd <- fd[is.finite(fd$uu) & fd$uu <100,]

## plot(fd)

## ## the beta vs coverage is shitty

## ## what if sampling them?
