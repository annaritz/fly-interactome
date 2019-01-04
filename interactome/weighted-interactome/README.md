# Weighting the Fly Interactome

This directory contains weighted interactome files of the form `Node1 Node2 Weight`.  This goes through the steps for the interactome collapsed by common name and merged evidence in FlyBase identifiers (`interactome-flybase-collapsed-evidence.txt`)


## Weighting Approach

We weight the edges using an equation with three terms.  The weighting is inspired by the [HIPPIE interactome](http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/index.php) weighting [described in Schaefer et al, Plos One, 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031826).  

1. _s1(studies)_: non-linear saturating function of the number of studies (PubMedIDs)
2. _s2(databases)_: non-linear saturating function of the number of databases that report interaction
3. _s3(evidence)_: Bayesian weighting function from [Yeger-Lotem et al., Nature Genetics 2009](https://doi.org/10.1038/ng.337) and also described in [Poirel et al., J. Computational Biology 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3646337/).

The final score is a weighted linear combination of the three terms, where the weights sum to 1:

_w1_ _s1(studies)_ + _w2_ _s2(databases)_ + _w3_ _s3(evidence)_

There are five unknowns in this equation: the steepness of the slope in the two non-linear saturating functions (call them _a1_ and _a2_), and _w1_ and _w2_ (since _w3 = 1-w2-w3_).  These are chosen via a parameter sweep using the DroID dataset (TODO).  This is similar to the process described in [described in Schaefer et al, Plos One, 2012](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031826).

### Evidence-based Weighting: 

#### Input Files

The file `weight-edges-by-evidence.py` is adapted from [Poirel et al., J. Computational Biology 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3646337/) to weight edges by GO term similarity.  It takes a number of required files:

1. The interactome (collapsed version of file default): `../interactome-flybase-collapsed-evidence.txt`

2. Gene Annotations file. Species-specific file from [Gene Ontology](http://www.geneontology.org/page/download-go-annotations).  See more information about the [FlyBase README](http://geneontology.org/gene-associations/readme/fb.README).  

```
wget http://geneontology.org/gene-associations/gene_association.fb.gz
gunzip gene_association.fb.gz
```

Convert the FlyBase-formatted gene association file to a `.gmt` file (this is the format that the script expects):

```
python ../../utils/generate_annotation_file.py gene_association.fb gene_association.gmt
```

3. Gene Ontology file (.obo).  The daily releases of the downloads (downloaded 2018-12-30):

```
wget http://purl.obolibrary.org/obo/go.obo
```

4. GO Functions for weighting. These should be GO terms that include regulatory (or signaling) interactions that will be considered positives. We will make this very general, with the goal to identify protein regulation:

```
GO:0005515	protein binding (note: not in go.obo -- this is ignored).
GO:0023051	regulation of signaling
GO:0009966	regulation of signal transduction (note: descendant of "regulation of signaling" -- this is ignored)
```

This is stored in `functions.txt`.

#### Running the `weight-edges-by-evidence.py` script

```
python weight-edges-by-evidence.py -n ../interactome-flybase-collapsed-evidence.txt -a gene_association.gmt -t go.obo -f functions.txt -o interactome-flybase-collapsed-evidence-weighted --maxsetsize 1000
```

To plot histograms of edges:
```
python ../../utils/plot-hist.py interactome-flybase-collapsed-evidence-weighted.txt 2 50 "Weighted Edges" "Edge Weight" "# of Edges" interactome-flybase-collapsed-evidence-weighted.png
python ../../utils/plot-hist.py interactome-flybase-collapsed-evidence-weighted-edge_type_weights.txt 3 20 "Weighted Edge Types" "Edge Type Weight" "# of Edge Types" interactome-flybase-collapsed-evidence-weighted-edge_type_weights.png
```

Once you have calculated probabilities for the entire network, you can quickly reweight subsets (which will be useful for the assessment section) by passing in one of the intermediate files (`*edge_type_probs.txt`):

```
python weight-edges-by-evidence.py -n ../interactome-flybase-collapsed-evidence.txt --probs interactome-flybase-collapsed-evidence-weighted-edge_type_probs.txt -o test_fast_run
```

### Weighting Edges

Usage:
```
Usage: python weight-edges.py [OPTIONS]

Options:
  -h, --help            show this help message and exit
  -n STR, --network=STR
                        protein-protein interactome (tab-delimited). Required.
  -c, --collapsed       Interactome is collapsed by common name (has 7
                        columns).
  -o STR, --outprefix=STR
                        output file prefix. Required.
  -e STR, --evidence=STR
                        Evidence weights (from weight-edges-by-evidence.py)
  --a1=FLOAT            Parameter that controls the steepness of s1 (study-
                        based term). Default = 1.
  --a2=FLOAT            Parameter that controls the steepness of s2 (db-based
                        term). Default = 1
  --w1=FLOAT            Weight of s1 (study-based term).  Float between 0 and
                        1; w3 = 1-w1-w2. Default = 0.33
  --w2=FLOAT            Weight of s2 (db-based term).  Float between 0 and 1;
                        w3 = 1-w1-w2. Default = 0.33
```

Test run:

```
python weight-edges.py -n ../interactome-flybase-collapsed-evidence.txt -c -e interactome-flybase-collapsed-evidence-weighted.txt -o weighted-test.txt
python ../../utils/plot-hist.py weighted-test.txt 2 50 "Final Edges (default params)" "Final Edge Weights" "# of Edges" weighted-test.png
```

### Parameter Sweep

```
python parameter-sweep.py ../interactome-flybase-collapsed-evidence.txt interactome-flybase-collapsed-evidence-weighted.txt interactome-flybase-collapsed-evidence-weighted-edge_type_probs.txt droid force
```

The last argument overwrites all files in the directories -- if it is missing, then files will only be generated if they are missing.


### NMII Details

Non-muscle myosin (NMII) is composed of mutliple subunits, each with their own protein name in FlyBase and [UniProt](https://www.uniprot.org/uniprot/?query=%22non%20muscle%22%20myosin&fil=organism%3A%22Drosophila+melanogaster+%28Fruit+fly%29+%5B7227%5D%22&sort=score):

- **zip**: zipper; Myosin heavy chain; CG15792 (UniProtID: Q99323)
- **Mlc-c**: Myosin essential light chain (UniProtID: P54357)
- **sqh**: Myosin regulatory light chain sqh; CG3595 (UniProtID: P40423)

I *think* that the regulatory regions are on the light chains -- Mlc-c (which binds myosin and calcium) and sqh (which binds calcium).  Let's start with Mlc-c -- it is annotated to 10 GO terms:

- [**GO:0016459**](http://amigo.geneontology.org/amigo/term/GO:0016459): myosin complex
- [**GO:0005509**](http://amigo.geneontology.org/amigo/term/GO:0005509): calcium ion binding
- [**GO:0032036**](http://amigo.geneontology.org/amigo/term/GO:0032036): myosin heavy chain binding
- [**GO:0030048**](http://amigo.geneontology.org/amigo/term/GO:0030048): actin filament-based movement
- [**GO:0031477**](http://amigo.geneontology.org/amigo/term/GO:0031477): myosin VII complex
- [**GO:0031476**](http://amigo.geneontology.org/amigo/term/GO:0031476): myosin VI complex
- [**GO:0031475**](http://amigo.geneontology.org/amigo/term/GO:0031475): myosin V complex
- [**GO:0017022**](http://amigo.geneontology.org/amigo/term/GO:0017022): myosin binding
- [**GO:0005515**](http://amigo.geneontology.org/amigo/term/GO:0005515): protein binding
- [**GO:0016460**](http://amigo.geneontology.org/amigo/term/GO:0016460): myosin II complex

### Developer Notes

The weighting code produces **different** weights than the original evidence-based weighting. The issue comes down to what's specified as the space of possible pairs -- the previous method took all (_n_ choose 2) pairs for the _n_ nodes in the network.  This is an issue because the majority of those pairs are not in the network; as a result, most of the sampled negatives were never considered because they weren't edges in the network.  Instead, this method defines the space as all _interacting_ pairs (all edges in the network).  
