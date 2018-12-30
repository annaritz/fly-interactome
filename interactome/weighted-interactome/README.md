# Weighting the Fly Interactome

This directory contains weighted interactome files of the form `Node1 Node2 Weight`.  


## Weighting Approach

We weight the edges using an equation with three terms.  The weighting is inspired by the [HIPPIE interactome](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/index.php) weighting (described in Schaefer et al, Plos One, 2012)[https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031826].  

1. _s1(studies)_: non-linear saturating function of the number of studies (PubMedIDs)
2. _s2(databases)_: non-linear saturating function of the number of databases that report interaction
3. _s3(evidence)_: Bayesian weighting function from [Yeger-Lotem et al., Nature Genetics 2009](https://doi.org/10.1038/ng.337) and also described in [Poirel et al., J. Computational Biology 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3646337/).

The final score is a weighted linear combination of the three terms, where the weights sum to 1:

_w1_ _s1(studies)_ + _w2_ _s2(databases)_ + _w3_ _s3(evidence)_

There are five unknowns in this equation: the sttepness of the slope in the two non-linear saturating functions (call them _a1_ and _a2_), and _w1_ and _w2_ (since _w3 = 1-w2-w3_).  These are chosen via a parameter sweep using the DroID dataset (TODO).
