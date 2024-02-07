# Fly Interactome

This fly interactome is compiled from a set of publicly-available sources (see the Databases section below). To cite this database, please see our publication which uses the unweighted, collapsed interactome:

Andy Zhao et al. [From network analysis to experimental validation: identification of regulators of non-muscle myosin II contractility using the folded-gastrulation signaling pathway](https://bmcmolcellbiol.biomedcentral.com/articles/10.1186/s12860-023-00492-3). _BMC Molecular and Cell Biology_ 2023 (volume 24 number 32).

##  Most Recent Interactome (Nov 2019)

- `interactome-flybase-collapsed-weighted.txt` ([link to file](interactome/weighted-interactome/interactome-flybase-collapsed-weighted.txt))

## Interactome History

The files include the following:

- `interactome-uniprot.txt` ([link to file](interactome/interactome-uniprot.txt)): Interactome in [UniProtKB IDs](http://www.uniprot.org/)
- `interactome-flybase.txt` ([link to file](interactome/interactome-flybase.txt)): Interactome in [FlyBase IDs](http://flybase.org/).

These files and code are released under the GPL-3.0 license.  Summary statistics about these interactomes is in the [interactome/](interactome/) directory.

**New as of July 2018:** collapsed versions of the interactions in common name are now available.  See details in the [interactome/](interactome/) directory.

- `interactome-uniprot-collapsed.txt` ([link to file](interactome/interactome-uniprot-collapsed.txt)): Interactome in [UniProtKB IDs](http://www.uniprot.org/)
- `interactome-flybase-collapsed.txt` ([link to file](interactome/interactome-flybase-collapsed.txt)): Interactome in [FlyBase IDs](http://flybase.org/).

**New as of December 2018:** interactomes with evidence codes mapped to [Molecular Interactions Controlled Vocabulary](https://www.ebi.ac.uk/ols/ontologies/mi) terms are now available.  See details in the [interactome/](interactome/) directory.

- `interactome-uniprot-collapsed-evidence.txt` ([link to file](interactome/interactome-uniprot-collapsed-evidence.txt)): Interactome in [UniProtKB IDs](http://www.uniprot.org/)
- `interactome-flybase-collapsed-evidence.txt` ([link to file](interactome/interactome-flybase-collapsed-evidence.txt)): Interactome in [FlyBase IDs](http://flybase.org/).

**New as of January 2019:** the FlyBase collapsed interactome has been weighted according to Schaefer et al.  See the details in the README file in the [weighted directory](interactome/weighted-interactome/). The file is listed as `final-weighted-interactome.txt`, and was used in PageRank (in the [algorithms directory](algorithms/)).  As of **November 2019**, the collapsed file with weighted edges is now available:

- `interactome-flybase-collapsed-weighted.txt` ([link to file](interactome/weighted-interactome/interactome-flybase-collapsed-weighted.txt))

## Databases
The interactome is built by combining these databases. FlyBase and UniProt interactomes exist for each database, which is parsed before the interactome is compiled.  Links point to the appropriate directory.
* [DroID](databases/DroID)
* [SignaLink2.0](databases/SignaLink)
* [mentha](databases/Mentha)
* [flyReactome](databases/flyReactome)
* [FlyMine](databases/flyMine)
* [myProteinNet](databases/myProteinNet)

## Dependencies for Generating the Interactome

The interactomes are provided in this repo; however, you may wish to re-generate some of the intermediate files, or re-build the interactome from scratch.  If this is the case, each database includes instructions as well as the quick start details below.

- Tested with Python version 2.7.10
- If you want to build the interactome from scratch, you need the following packages:
 - [networkx](https://networkx.org/documentation/networkx-1.11/) version 1.11.
 - [mygene](https://pypi.python.org/pypi/mygene) Python client that interfaces with [MyGene.info](http://mygene.info/).  This is required for mapping IDs.
 - [intermine](http://intermine.readthedocs.org/en/latest/web-services/) is a open source data warehouse build specifically for the integration and analysis of complex biological data.  This is required to query FlyMine interactions (though these can also be acquired through the web interface).

For instructions for running the full interactome, see the instructions in the [interactomes/](interactomes/) directory.

## Developer Notes and TODOs
- Ported over from bitbucket account at research/2016-06-17-apical-constriction/interactomes
- [Flynet](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2773252/) is a transcriptional regulatory network.
- Add the [protein complex network](http://www.sciencedirect.com/science/article/pii/S0092867411010804), conducted in S2R+ cells.
- Add a requirements.yaml file for conda environment. Include networkx, matplotlib v2.2.
