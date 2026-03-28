# Between Cell Differentiation

This section explores how well Shannon entropy captures heterogeneity across features generated with varying levels of diversity.
We begin by generating a fixed number of Markov chains using a two-state transition matrix of the form:

$$
\begin{bmatrix}
\rho & 1 - \rho \\
\rho & 1 - \rho
\end{bmatrix}
$$

with different values of $\rho$.

Each feature is then simulated by sampling from a subset of these Markov chains.
The size of this subset defines the heterogeneity score of the feature: features with higher heterogeneity scores are sampled from a larger pool of chains, and therefore tend to exhibit higher Shannon entropy.
Conversely, features with lower heterogeneity scores are drawn from fewer chains, resulting in lower entropy.
