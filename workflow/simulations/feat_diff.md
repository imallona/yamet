# Feature Differentiation

This section focuses on the ability of sample entropy to distinguish between features generated from markov chains with different transition matrices. We have 4 types of features and we sample `N` features, each of length `f` to create a template. We then sample from this template multiple times to simulate cell samples. Each feature has a starting distribution for the markov chain that is shared across samples.
