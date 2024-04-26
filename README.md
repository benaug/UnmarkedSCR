# UnmarkedSCR
Various unmarked SCR MCMC samplers in nimble

Currently, 2D (individual by trap) versions of the Poisson, Negative Binomial, and Bernoulli observation model are provided. There is also a version using presence/absence data (Ramsey).

4/26/24 Disclaimer: I found an error in the proposal probabilities for the version using Bernoulli data. I thought you could enforce the constraints stemming from the fact that each individual can only be detected once per trap-occasion using the 2D individual by trap data, but you must use the 3D data. I now use the 2D data in the nimble model file for all updates except y/ID, reconstruct the 3D data and use that to update y/ID, then convert it back to 2D to feed back to nimble. This is faster than using the 3D data in nimble, but does not allow occasion effects on detection. This error evaded my notice because the estimates were good in the low home range overlap scenarios I was using to test the code.