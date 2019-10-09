#!/usr/bin/env R

WD = '/home/imallona/src/meth_entropies/runs/tuple_lengths/EGAR00001733321_ind_43_N656_15_fv'

## only chr19
fns <- list.files(WD, '.*entropy.*19.bed')

d <- list()

for (fn in fns) {
    d[[fn]] <- read.table(file.path(WD, fn),
                          header = FALSE)
    colnames(d[[fn]]) <- c('chr', 'start', 'end', 'id', 'entropy', 'strand', 'methylation')
}


colors <- c('black', 'forestgreen', 'brown', 'blue')

png('a2.png', width = 2000, height = 300)

plot(x = d[[1]]$start, y = d[[1]]$entropy, col = colors[1])
for (i in 2:length(d)) {
    lines(x = d[[i]]$start, y = d[[i]]$entropy, col = colors[i], type = 'p')
}

dev.off()


png('a3.png', width = 2000, height = 300)

plot(x = head(d[[1]]$start, 1000), y = head(d[[1]]$entropy, 1000),
     pch = 19,
     col = colors[1], type = 'p')
for (i in 2:length(d)) {
    lines(x = head(d[[i]]$start, 1000), y = head(d[[i]]$entropy, 1000),
          pch = 19,
          col = colors[i], type = 'p')
}

dev.off()


par(mfrow = c(4,1))

for (i in 1:length(d)) {
    plot(x = head(d[[i]]$start, 1000), y = head(d[[i]]$entropy, 1000),
         pch = 19,
         col = colors[i], type = 'p',
         xlab = names(d)[i])
}

