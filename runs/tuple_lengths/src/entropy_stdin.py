#!/bin/python3
##
## Calculates entropies and methylations per row for whathever the output
##   of methtuple it is (regardless of the -m value)
##
## Reads from stdin, writes to stdout
## @todo optimize
##
## Izaskun Mallona
## 08 Oct 2019

import sys
from scipy.stats import entropy
import pandas

# read from stdin as many columns as methtuple readouts

fd = pandas.read_csv(sys.stdin, skiprows = 0, sep = '\t', header = 0)

fd2 = fd.filter(regex=r'^[MU].*', axis=1)
starter = fd.filter(regex=r'^pos.*', axis=1).head()

# calculate probabilities for each item (not needed)
# fd2 = fd.div(fd.sum(axis=1), axis=0)

# calculate entropy
entropies = fd2.apply(entropy, axis=1)

# get the proportion of `M` elements within each tuple, to leverage the contribution
# to methylation
tmp = fd2.columns.str.count('M')/len(fd2.columns.values[0])
fd3 = fd2.multiply(tmp, axis = 'columns')      
methylations = fd3.apply(sum, axis = 1)/fd2.apply(sum, axis = 1)

# pasting everything back: bedfile
# originally from awk snippet (2-tuples)
# print $1,$3,$4,"MM"$5";MU"$6";UM"$7";UU"$8,E,$2

## currently, the tuple has no name
empty_name = pandas.Series()

print(pandas.concat([fd[['chr']],
                     fd[['pos1']],
                     fd[starter.columns[-1:]],
                     empty_name,
                     entropies,
                     fd[['strand']],
                     entropies, methylations],
                    axis=1,
                    ignore_index = True).to_string(index = False, header = False,
                                                   justify = 'unset'))

# df.to_csv('dfsavename.csv.gz', compression='gzip')
