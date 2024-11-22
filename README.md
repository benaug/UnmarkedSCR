# UnmarkedSCR
Various unmarked SCR MCMC samplers in nimble

There are two data augmentation approaches, 1) regular data augmentation as used in the original paper and 2) an alternative that allows a Poisson prior on expected abundance (DA2 in file names). The latter allows faster N/z updates. Multisession versions use DA2.

There are 5 observation models, 1) Poisson, 2) Bernoulli, 3) negative binomial, 4) zero-truncated Poisson hurdle, and 5) Bernoulli-Ramsey. Currently, 3 and 5 are only available with regular data augmentation.

See test scripts and then the files loaded inside test scrips, e.g., model files, data simulator, data initializer, nimble functions/custom updates. Generally, these models are garbage, but they serve as a basis for better models that use partial ID information.

Perhaps a good place to start is to check out test scripts and model files for:
Poisson DA2, Bernoulli DA2, hurdleZTPois DA2. These are set up for the same scenario with an informative prior on sigma.


4/26/24 Disclaimer: I found an error in the proposal probabilities for the version using Bernoulli data. I thought you could enforce the constraints stemming from the fact that each individual can only be detected once per trap-occasion using the 2D individual by trap data, but you must use the 3D data. I now use the 2D data in the nimble model file for all updates except y/ID, reconstruct the 3D data and use that to update y/ID, then convert it back to 2D to feed back to nimble. This is faster than using the 3D data in nimble, but does not allow occasion effects on detection. This error evaded my notice because the estimates were good in the low home range overlap scenarios I was using to test the code.

11/21/24: I added a version of Poisson DA2 that uses the observation model likelihood that is marginalized over individuals.
Further, I use results from Herliansyah et al. (2024) to speed up N/z and s updates.
https://link.springer.com/article/10.1007/s13253-023-00598-3


See testscript for Poisson DA2 Marginal or Poisson DA2 Marginal Multisession