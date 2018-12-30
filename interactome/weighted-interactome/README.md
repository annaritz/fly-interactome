# Weighting the Fly Interactome

This directory contains weighted interactome files of the form `Node1 Node2 Weight`.  


## Weighting Approach

We weight the edges using an equation with three terms.  The weighting is inspired by the [HIPPIE interactome](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/index.php) weighting (described in Schaefer et al, Plos One, 2012)[https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031826].  

1. $s_1$: non-linear saturating function of the number of studies (PubMedIDs)
2. Term 2: non-linear saturating function of the number of databases that report interaction
3. Term 3: Bayesian weighting function from [Yeger-Lotem et al., Nature Genetics 2009](https://doi.org/10.1038/ng.337) and also described in [Poirel et al., J. Computational Biology 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3646337/).

