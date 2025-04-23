# Within Cell Differentiation

This section focuses on a simulation setup that isolates within-cell sequence entropy from across-cell entropy.
We begin with two sequence templates and generate an equal number of copies of each, arranged sequentially by group.
Each feature is assigned a heterogeneity score, which governs the degree of shuffling - highly heterogeneous features undergo greater shuffling of copies within the feature, while less heterogeneous features experience minimal rearrangement.
This controlled shuffling preserves the shannon entropy but modulates the sequence regularity.
Finally, we randomly select 30% of the positions from every feature of every cell to add heterogeneity to the average methylation level which is important for the comparison with scMET.
