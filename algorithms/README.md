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
