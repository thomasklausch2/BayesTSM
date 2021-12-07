# bayesTSM

**Bayesian 3-state AFT model for screening outcomes.**

Please note:

* This 'package' currently is a collection of functions. We plan to publish an `R` `library` soon.
* The main function is `bts_survreg`
* Use Rstudio project or define working directory path
* To load all functions required  `source` path `/controllers/startup.r`
* For a demonstration of the Bayesian model estimation by Gibbs check the `demo_MCMC.r` with demonstrations of `bts_survreg` using simulated data
* For a demonstration of the EM algorithm for the lognormal-lognormal model, see `demo_EM.r`
* The `foreach` package is used for parallel and sequential processing. If your system is incompatible with `foreach` consider using the `bts_survreg_seq` function instead which uses a `for` loop insead (can be slow).

**Author information:** Thomas Klausch (t.klausch@amsterdamumc.nl). Reference as: Klausch, T., Akwiwu, E. U., van de Wiel, M. A., Coupé, V. M. H., Berkhof, J. (2021). A Bayesian accelerated failure time model for interval censored three-state screening outcomes. Available under: https://arxiv.org/abs/2110.02649

**License**: This product is provided under a creative commons license: Attribution-NonCommercial-ShareAlike
CC BY-NC-SA.

**Important note**: Throughout the code and output, the annotation `X` is used for the time from state 1 to state 2 and `S` is used for the time from state 2 to state 3. In the accompanying main paper the notation $x$ and $t$ are used for these times, respectively. Letters `t` and `T` are problematic to use in `R` as they have frequently used other definitions. To avoid conflicts the notation `S` is used throughout.
