# Algorithms

These algorithms use Python's [Networkx v1.11](https://networkx.github.io/documentation/networkx-1.11/) (TODO - upgrade to 2.x).  

## Teleporting Random Walk (PageRank)

Computes a random walk on an undirected graph (?) from a set of source nodes, teleporting back to one of the source nodes uniformly at random with probability _q_.  This implementation of PageRank is modified from "Top-down network analysis to drive bottom-up modeling of physiological processes" ([Poirel et al., J. Comp. Biol, 2013]((https://www.ncbi.nlm.nih.gov/pubmed/23641868))).

For example, to compute the random walk from the essential light chain of NMII (Mlc-c) with different teleportation probabilities:

```
python PageRank.py --undirected --output pagerank_mlc-c_q_0.2.txt -q 0.2 -t Mlc-c ../interactome/weighted-interactome/final-weighted-interactome.txt
python PageRank.py --undirected --output pagerank_mlc-c_q_0.5.txt -q 0.5 -t Mlc-c ../interactome/weighted-interactome/final-weighted-interactome.txt
python PageRank.py --undirected --output pagerank_mlc-c_q_0.8.txt -q 0.8 -t Mlc-c ../interactome/weighted-interactome/final-weighted-interactome.txt
``` 

## Visualize network.

This posts the top 100 nodes in the PageRank output to GraphSpace, and annotates the nodes and edges.  Warning: there are two hard-coded annotation files required, and GraphSpace module is also required. Really need to make a virtual environment for this...

```
python3 ../utils/viz-network.py ../interactome/interactome-flybase-collapsed-evidence.txt ../interactome/weighted-interactome/final-weighted-interactome.txt pagerank_mlc-c_q_0.5.txt 0.5
```

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
