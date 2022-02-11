# UnmarkedSCR
Various unmarked SCR MCMC samplers in nimble

Currently, 2D (individual by trap) versions of the Poisson and Bernoulli observation model are provided. The latter is probably rarely appropriate, and is generally even worse in terms of precision and bias than the Poisson version (count data is just more informative than binary data). The binary version is also slower because when using the 2D data, we must check that we never assign 2 samples from the same trap-occasion to the same individual during the ID update.