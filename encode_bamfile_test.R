#!/usr/bin/env R
##
## Cometh agreement scoring simulations
##
## Izaskun Mallona
## 12th nov 2018

library('latex2exp')
library(entropy)
library(Cairo)
## library(MASS)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

NROWS <- 1e6
MINCOVERAGE <- 5 ## this is enforced at the bash script!

WD <- file.path('/home/imallona', 'mnt', 'nfs', 'cg_shadows', 'data', 'ENCFF857QML')

ssample <- 'ENCFF857QML'

meth_fn <- file.path(WD, sprintf('%s.CG.2_meth_colored.bed', ssample))
entropy_fn <- file.path(WD, sprintf('%s.CG.2_entropy_colored.bed', ssample))

beta2m <- function(beta) {
    m <- log2(beta/(1 - beta))
    m
}

add.alpha <- function(col, alpha=1){
 
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) rgb(x[1], x[2], x[3], alpha = alpha))  
}


d <- list (entropy = read.table(entropy_fn,
                                colClasses = c("NULL", "integer", "integer", "factor", "numeric",
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


                                
colnames(d$entropy) <- c('start', 'end', 'code', 'entropy', 'hmm')
colnames(d$meth) <- c('beta', 'hmm')
d$meth$beta <- d$meth$beta/1000
d$entropy$distance <- d$entropy$end - d$entropy$start
## d$entropy$mm <- sapply(strsplit(as.character(d$entropy$code), ';'),
##                        function(x) return(gsub('MM', '', x)))


tmp  <- sapply(strsplit(as.character(d$entropy$code), ';'),
               function(x) return(substr(x,3, nchar(x))))

d$entropy$mm <- as.numeric(tmp[1,])
d$entropy$mu <- as.numeric(tmp[2,])
d$entropy$um <- as.numeric(tmp[3,])
d$entropy$uu <- as.numeric(tmp[4,])

rm(tmp)

d$entropy$coverage <- d$entropy$mm + d$entropy$um + d$entropy$mu + d$entropy$uu

## d <- d[d$coverage > MINCOVERAGE,]

for (color in colors) {
    d$entropy[,color] <- NA
    d$entropy[,color] <- grepl(color, as.character(d$entropy$hmm))
}



d$entropy$in_boundary <- apply(d$entropy[,as.character(colors)], 1, function(x) sum(x) > 1)
table(d$entropy$in_boundary)


d$entropy$beta <- d$meth$beta 
d <- d$entropy
d$log10dist <- log10(d$distance)


targets <- c('coverage',
             'distance',
             'beta',
             'entropy')


## test standardization start

## test standardization end


png(file.path(WD, sprintf('%s_exploratory_1.png', ssample)), width = 700,
         height = 700)

par(cex.axis = 1.4,
    cex.lab = 1.4,
    cex.main = 1.4,
    cex.sub = 1.4,
    pty = "s",
    mar=c(5.1,4.1,4.1,2.1),
    oma = c(1, 0, 0, 0))

plot(d[,targets], pch = 20,
     ## col = add.alpha(ifelse(d$in_boundary, 'darkred', 'black' ), 0.5))
     col = add.alpha('black', 0.5))
dev.off()


targets <-  c('coverage',
             'log10dist',
             'beta',
             'entropy')

entropy_max <- max(d$entropy)


rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)


for (color in colors) {

    curr <- d[d[,color],]
    png(file.path(WD, sprintf('exploratory_1_%s.png', color)),
        width = 700,
        height = 700)

    plot(curr[,targets],
         pch = 20,
         col = add.alpha(ifelse(curr[,"in_boundary"], 'red', 'black'), 0.5),
         main = color)
    
    dev.off()

    
    png(file.path(WD, sprintf('exploratory_2_%s.png', color)),
        width = 400,
        height = 400)
        
    plot(x = curr$beta,
         y = curr$entropy,
         xlim = c(0,1),
         ylim = c(0, entropy_max),
         pch = 20,
         col = add.alpha(ifelse(curr[,"in_boundary"], 'red', 'black'), 0.5),
         main = color)
    
    dev.off()

    p <- ggplot(curr[,c('beta', 'entropy')], aes(beta, entropy))

    ## h3 <- p + stat_bin2d(bins=25) + scale_fill_gradientn(colours=r) + xlim(-0.1, 1.1) +
    ##     ylim(-0.1, entropy_max + 0.1) + ggtitle(color)
    
    h3 <- p + stat_bin2d(bins=25) + scale_fill_gradientn(colours=r) + xlim(-0.1, 1.1) +
        ylim(-0.1, entropy_max + 0.1) + ggtitle(color) + theme(text = element_text(size=20)) +
        xlab('DNA methylation (beta)') + ylab("Shannon's entropy (H)")
    
    
    ggsave(h3, file = file.path(WD, sprintf('exploratory_3_%s.png', color)))
    
  
}


## here ################################################################################

## stacked boxplots

## which is the proportion of 0 entropy for each HMM discarding boundaries?



d$unmethylated <- d$beta < 0.2
d$methylated <- d$beta >= 0.8
d$zero_entropy <- d$entropy == 0


non_boundary <- list(median = tapply(d[!d$in_boundary, 'entropy'],
                                     as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                     function(x) quantile(x, probs = 0.5)),
                     lower = tapply(d[!d$in_boundary, 'entropy'],
                                     as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                  function(x) quantile(x, probs = 0.25)),
                     upper = tapply(d[!d$in_boundary, 'entropy'],
                                     as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                  function(x) quantile(x, probs = 0.755)))

## sort by median, plot
sorted <- names(sort(non_boundary$median))

png('quantiles_entropy_no_boundaries.png', height = 480, width = 550)

par(cex.axis = 1.8,
    cex.lab = 1.8,
    cex.main = 1.8,
    cex.sub = 1.8,
    pty = "s",
    mar=c(5.1,4.1,4.1,2.1),
    oma = c(1.2, 0, 0, 0),
    xpd = TRUE)

plot(non_boundary$lower[sorted], pch = 20,
     ylim = c(0, max(unlist(non_boundary))),
     ylab = "Shannon's entropy (H)" ,
     xaxt = "n",
     xlab = "",
     main = '')

points(non_boundary$median[sorted], pch = 20, col = 'red')
points(non_boundary$upper[sorted], pch = 20)

segments(y0 = non_boundary$lower[sorted],
         y1 = non_boundary$upper[sorted],
         x0 = 1:length(non_boundary$lower),
         x1 = 1:length(non_boundary$lower))

legend('topright', col = c('black', 'red', 'black'), pch = 20, c('3Q', 'median', '1Q'),
       inset=c(-0.25,0))

axis(1, at = 1:length(sorted), labels = sorted, las = 2)

dev.off()



non_boundary <- list(median = tapply(d[!d$in_boundary, 'beta'],
                                     as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                     function(x) quantile(x, probs = 0.5)),
                     lower = tapply(d[!d$in_boundary, 'beta'],
                                     as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                  function(x) quantile(x, probs = 0.25)),
                     upper = tapply(d[!d$in_boundary, 'beta'],
                                     as.factor(as.character(d[!d$in_boundary, 'hmm'])),
                                  function(x) quantile(x, probs = 0.755)))



png('quantiles_methylation_no_boundaries.png')

par(cex.axis = 1.8,
    cex.lab = 1.8,
    cex.main = 1.8,
    cex.sub = 1.8,
    pty = "s",
    mar=c(5.1,4.1,4.1,2.1),
    oma = c(1.2, 0, 0, 0))

plot(non_boundary$lower[sorted], pch = 20,
     ylim = c(0, max(unlist(non_boundary))),
     ylab = 'methylation (beta)',
     xaxt = "n",
     xlab = "",
     main = '')

points(non_boundary$median[sorted], pch = 20, col = 'red')
points(non_boundary$upper[sorted], pch = 20)
segments(y0 = non_boundary$lower[sorted],
         y1 = non_boundary$upper[sorted],
         x0 = 1:length(non_boundary$lower),
         x1 = 1:length(non_boundary$lower))

axis(1, at = 1:length(sorted), labels = sorted, las = 2)

dev.off()



## stratified by dnameth (low, mid, high, cutoffs 0.2 and 0.8 beta)


curr <- d[!d$in_boundary,]
curr$hmm <- as.factor(as.character(curr$hmm))
## categorical methylation status
curr$meth_cat <- 'mid'
curr$meth_cat[curr$unmethylated] <- 'low'
curr$meth_cat[curr$methylated] <- 'high'
curr$meth_cat <- factor(curr$meth_cat, levels = c('low', 'mid', 'high'))

## ggplot(curr, aes(factor(hmm), entropy, fill = meth_cat)) + 
##   geom_bar(stat="identity", position = "dodge") + 
##   scale_fill_brewer(palette = "Set1")



h <- ggplot(aes(y = entropy, x = hmm, fill = meth_cat), data = curr) +
    geom_boxplot(outlier.alpha = 0.5, color = 'darkblue') +
    scale_fill_manual(values=c("gray90", "gray60", "gray30")) +
    ylab("Shannon's entropy (H)") +
    xlab('chromatin state (chromHMM)') +
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(h, file = file.path('entropy_no_boundary_boxplot_stratified_by_meth.png'),
       height = 4, width = 8)


## how many instances are there?

## ggplot(curr, aes(factor(hmm),fill = meth_cat)) +
##     geom_bar(stat="count", position = "dodge") + 
##     scale_fill_manual(values=c("gray90", "gray60", "gray30")) +
##     theme_bw() + 
##     theme(text = element_text(size = 15),
##           axis.text.x = element_text(angle = 90, hjust = 1))


## h1 <- ggplot(curr, aes(factor(hmm),fill = meth_cat)) +
##     scale_y_log10(
##         breaks = scales::trans_breaks("log10", function(x) 10^x),
##         labels = scales::trans_format("log10", scales::math_format(10^.x))
##     ) +    
##     geom_bar(stat="count", position = "dodge", color = 'black') + 
##     scale_fill_manual(values=c("gray90", "gray60", "gray30")) +
##     theme_bw() + 
##     theme(text = element_text(size = 15),
##           axis.text.x = element_text(angle = 90, hjust = 1)) +
##     annotation_logticks(sides = 'l') 

h1 <- ggplot(curr, aes(factor(hmm),fill = meth_cat)) +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +    
    geom_bar(stat="count", width = 0.6, position = position_dodge(width=0.7), color = 'black') + 
    scale_fill_manual(values=c("gray90", "gray60", "gray30")) +
    theme_bw() +
    xlab('chromatin state (chromHMM)') +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    annotation_logticks(sides = 'l') 


ggsave(h1, file = file.path('counts_no_boundary_boxplot_stratified_by_meth.png'),
       height = 4, width = 8)


## same with those in boundaries!
## (but grouping by single colors, not the pairwise relationships)

curr <- d[d$in_boundary,]
fd <- melt(curr[c('entropy', 'beta', as.character(colors))],
           c("entropy", 'beta'),
           variable.name = "hmm")
fd <- fd[fd$value,]


## categorical methylation status
fd$meth_cat <- 'mid'
fd$meth_cat[fd$beta < 0.2] <- 'low'
fd$meth_cat[fd$beta >= 0.8] <- 'high'
fd$meth_cat <- factor(fd$meth_cat, levels = c('low', 'mid', 'high'))



h <- ggplot(aes(y = entropy, x = hmm, fill = meth_cat), data = fd) +
    geom_boxplot(outlier.alpha = 0.5) +
    scale_fill_manual(values=c("gray90", "gray60", "gray30")) +
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(h, file = file.path('entropy_boundaries_boxplot_stratified_by_meth.png'))


h1 <- ggplot(fd, aes(factor(hmm),fill = meth_cat)) +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +    
    geom_bar(stat="count", width = 0.6, position = position_dodge(width=0.7), color = 'black') + 
    scale_fill_manual(values=c("gray90", "gray60", "gray30")) +
    theme_bw() + 
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    annotation_logticks(sides = 'l') 


ggsave(h1, file = file.path('counts_boundaries_boxplot_stratified_by_meth.png'))



## on some reldists to marks
## without cpgislands

reld <- list()

for (fn in list.files(file.path(WD, 'reldist_no_cpgi'), pattern = "*reldist$")) {
    reld[[fn]] <- read.table(file.path(WD, 'reldist_no_cpgi', fn),
                    sep = "\t", header = TRUE)
}




png( sprintf('dnasei_narrow_reldistplot.png', fn),
    width = 480, height = 480)
par(mfrow = c(1,1),
    oma = c(5,4,1,1) + 0.1,
    mar = c(2,1,2,1) + 0.1,
    pty = "s",
    cex.main = 1.4,
    cex.lab = 1.4,
    cex.axis = 1.4)

## d <- read.table(fn,
##                 sep = "\t", header = TRUE)

plot(reld$ENCFF857QML_lowly_meth_low_entropy__vs__ENCFF030XPN_dnasei_narrow.bed.reldist$reldist,
     reld$ENCFF857QML_lowly_meth_low_entropy__vs__ENCFF030XPN_dnasei_narrow.bed.reldist$fraction,
     type = 'b',
     pch = 1,
     ## main = sprintf('%s vs\n %s', tag_a, tag_b),
     axes = TRUE,
     xlab = "",
     ylab = "",
     ylim = c(0, 1),
     col = 'black',
     main = 'Entropy vs DNAseI in hESC H1')

points(reld$ENCFF857QML_highly_meth_low_entropy__vs__ENCFF030XPN_dnasei_narrow.bed.reldist$reldist,
     reld$ENCFF857QML_highly_meth_low_entropy__vs__ENCFF030XPN_dnasei_narrow.bed.reldist$fraction,
     type = 'b',
     pch = 2,
     col = 'darkred')



points(reld$ENCFF857QML_lowly_meth_high_entropy__vs__ENCFF030XPN_dnasei_narrow.bed.reldist$reldist,
     reld$ENCFF857QML_lowly_meth_high_entropy__vs__ENCFF030XPN_dnasei_narrow.bed.reldist$fraction,
     type = 'b',
     pch = 3,
     col = 'darkgreen')

points(reld$ENCFF857QML_highly_meth_high_entropy__vs__ENCFF030XPN_dnasei_narrow.bed.reldist$reldist,
     reld$ENCFF857QML_highly_meth_high_entropy__vs__ENCFF030XPN_dnasei_narrow.bed.reldist$fraction,
     type = 'b',
     pch = 4,
     col = 'darkblue')


legend('topright',
       col = c('black', 'darkred', 'darkgreen', 'darkblue'),
       pch = 1:4,
       c('low meth low entropy', 'high meth low entropy', 'low meth high entropy', 'high meth high entropy '))

title(xlab = sprintf("\nrelative distance to DNAse I peaks ENCFF857QML\noutside CpG islands"),
      ylab = "frequency",
      outer = TRUE, line = 2)

dev.off()


## what about ctcf?


png( sprintf('ctcf_reldistplot.png', fn),
    width = 480, height = 480)
par(mfrow = c(1,1),
    oma = c(5,4,1,1) + 0.1,
    mar = c(2,1,2,1) + 0.1,
    pty = "s",
    cex.main = 1.4,
    cex.lab = 1.4,
    cex.axis = 1.4)

## d <- read.table(fn,
##                 sep = "\t", header = TRUE)

plot(reld$ENCFF857QML_lowly_meth_low_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$reldist,
     reld$ENCFF857QML_lowly_meth_low_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$fraction,
     type = 'b',
     pch = 1,
     ## main = sprintf('%s vs\n %s', tag_a, tag_b),
     axes = TRUE,
     xlab = "",
     ylab = "",
     ylim = c(0, 1),
     col = 'black',
     main = 'Entropy vs CTCF in hESC H1')

points(reld$ENCFF857QML_highly_meth_low_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$reldist,
     reld$ENCFF857QML_highly_meth_low_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$fraction,
     type = 'b',
     pch = 2,
     col = 'darkred')



points(reld$ENCFF857QML_lowly_meth_high_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$reldist,
     reld$ENCFF857QML_lowly_meth_high_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$fraction,
     type = 'b',
     pch = 3,
     col = 'darkgreen')

points(reld$ENCFF857QML_highly_meth_high_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$reldist,
     reld$ENCFF857QML_highly_meth_high_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$fraction,
     type = 'b',
     pch = 4,
     col = 'darkblue')


legend('topright',
       col = c('black', 'darkred', 'darkgreen', 'darkblue'),
       pch = 1:4,
       c('low meth low entropy', 'high meth low entropy', 'low meth high entropy', 'high meth high entropy '))

title(xlab = sprintf("\nrelative distance to CTCF peaks ENCFF368LWM\noutside CpG islands"),
      ylab = "frequency",
      outer = TRUE, line = 2)

dev.off()

## rampage

png( sprintf('rampage_reldistplot.png', fn),
    width = 480, height = 480)
par(mfrow = c(1,1),
    oma = c(5,4,1,1) + 0.1,
    mar = c(2,1,2,1) + 0.1,
    pty = "s",
    cex.main = 1.4,
    cex.lab = 1.4,
    cex.axis = 1.4)

## d <- read.table(fn,
##                 sep = "\t", header = TRUE)

plot(reld$ENCFF857QML_lowly_meth_low_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$reldist,
     reld$ENCFF857QML_lowly_meth_low_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$fraction,
     type = 'b',
     pch = 1,
     ## main = sprintf('%s vs\n %s', tag_a, tag_b),
     axes = TRUE,
     xlab = "",
     ylab = "",
     ylim = c(0, 1),
     col = 'black',
     main = 'Entropy vs Rampage in hESC H1')

points(reld$ENCFF857QML_highly_meth_low_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$reldist,
     reld$ENCFF857QML_highly_meth_low_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$fraction,
     type = 'b',
     pch = 2,
     col = 'darkred')



points(reld$ENCFF857QML_lowly_meth_high_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$reldist,
     reld$ENCFF857QML_lowly_meth_high_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$fraction,
     type = 'b',
     pch = 3,
     col = 'darkgreen')

points(reld$ENCFF857QML_highly_meth_high_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$reldist,
     reld$ENCFF857QML_highly_meth_high_entropy__vs__ENCFF368LWM_ctcf.bed.reldist$fraction,
     type = 'b',
     pch = 4,
     col = 'darkblue')


legend('topright',
       col = c('black', 'darkred', 'darkgreen', 'darkblue'),
       pch = 1:4,
       c('low meth low entropy', 'high meth low entropy', 'low meth high entropy', 'high meth high entropy '))

title(xlab = sprintf("\nrelative distance to Rampage peaks ENCFF788SVK\noutside CpG islands"),
      ylab = "frequency",
      outer = TRUE, line = 2)

dev.off()

## kdm1a

png( sprintf('kdm1a_reldistplot.png', fn),
    width = 480, height = 480)
par(mfrow = c(1,1),
    oma = c(5,4,1,1) + 0.1,
    mar = c(2,1,2,1) + 0.1,
    pty = "s",
    cex.main = 1.4,
    cex.lab = 1.4,
    cex.axis = 1.4)

## d <- read.table(fn,
##                 sep = "\t", header = TRUE)

plot(reld$ENCFF857QML_lowly_meth_low_entropy__vs__ENCFF562OAN_kdm1a.bed.reldist$reldist,
     reld$ENCFF857QML_lowly_meth_low_entropy__vs__ENCFF562OAN_kdm1a.bed.reldist$fraction,
     type = 'b',
     pch = 1,
     ## main = sprintf('%s vs\n %s', tag_a, tag_b),
     axes = TRUE,
     xlab = "",
     ylab = "",
     ylim = c(0, 1),
     col = 'black',
     main = 'Entropy vs KM1A in hESC H1')

points(reld$ENCFF857QML_highly_meth_low_entropy__vs__ENCFF562OAN_kdm1a.bed.reldist$reldist,
     reld$ENCFF857QML_highly_meth_low_entropy__vs__ENCFF562OAN_kdm1a.bed.reldist$fraction,
     type = 'b',
     pch = 2,
     col = 'darkred')



points(reld$ENCFF857QML_lowly_meth_high_entropy__vs__ENCFF562OAN_kdm1a.bed.reldist$reldist,
     reld$ENCFF857QML_lowly_meth_high_entropy__vs__ENCFF562OAN_kdm1a.bed.reldist$fraction,
     type = 'b',
     pch = 3,
     col = 'darkgreen')

points(reld$ENCFF857QML_highly_meth_high_entropy__vs__ENCFF562OAN_kdm1a.bed.reldist$reldist,
     reld$ENCFF857QML_highly_meth_high_entropy__vs__ENCFF562OAN_kdm1a.bed.reldist$fraction,
     type = 'b',
     pch = 4,
     col = 'darkblue')


legend('topright',
       col = c('black', 'darkred', 'darkgreen', 'darkblue'),
       pch = 1:4,
       c('low meth low entropy', 'high meth low entropy', 'low meth high entropy', 'high meth high entropy '))

title(xlab = sprintf("\nrelative distance to KM1A peaks ENCFF562OAN\noutside CpG islands"),
      ylab = "frequency",
      outer = TRUE, line = 2)

dev.off()


## without














## what about checking preferential boundaries, the ones that are the most frequent?

stop('till here')




## # Log scaling
## h3 <- p + stat_bin2d(bins=25) + scale_fill_gradientn(colours=r, trans="log")
## h3



## wide to long?

## summary(d$entropy[!d$entropy$in_boundary, 'entropy'])

## fit1 <- aov(entropy~hmm, data = d$entropy)
## summary(fit1)

stripchart(d$entropy$entropy~ d$entropy$hmm, jitter = 0.1,
           vertical = TRUE,
           pch = 21)


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
d$m <- beta2m(d$beta)

## targets <- c('coverage', 'distance', 'beta',
##              'm', 'score1', 'score2',
##              'entropy_four', 'entropy_three')

targets <- c('coverage',
             'distance', 'log10dist',
             'beta',
             'm',
             'entropy')

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
