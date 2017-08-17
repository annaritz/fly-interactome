# fly-interactome
Fly interactome compiled from a set of publicly-available sources.  The files include the following:

- `interactome-uniprot.txt` ([link to file](interactome/interactome-uniprot.txt)): Interactome in [UniProtKB IDs](http://www.uniprot.org/)
- `interactomeflybase.txt` ([link to file](interactome/interactome-flybase.txt)): Interactome in [FlyBase IDs](http://flybase.org/).

These files and code are released under the GPL-3.0 license.  Summary statistics about these interactomes is in the [interactome/](interactome/) directory.

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
 - [mygene](https://pypi.python.org/pypi/mygene) Python client that interfaces with [MyGene.info](http://mygene.info/).  This is required for mapping IDs.
 - [intermine](http://intermine.readthedocs.org/en/latest/web-services/) is a open source data warehouse build specifically for the integration and analysis of complex biological data.  This is required to query FlyMine interactions (though these can also be acquired through the web interface).

For instructions for running the full interactome, see the instructions in the [interactomes/](interactomes/) directory.

## Development Notes
- Ported over from bitbucket account at research/2016-06-17-apical-constriction/interactomes
- [Flynet](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2773252/) is a transcriptional regulatory network.