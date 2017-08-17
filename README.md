# fly-interactome
Fly interactome compiled from a set of publicly-available sources.

## Dependencies
- Tested with Python version 2.7.10
- If you want to build the interactome from scratch, you need the [mygene](https://pypi.python.org/pypi/mygene) Python client that interfaces with [MyGene.info](http://mygene.info/).  This is required for mapping IDs.

## Quick Start

This code outputs two interactomes: one interactome based on [UniProtKB identifiers](http://www.uniprot.org/), and one interactome base on [FlyBase identifiers](http://flybase.org/).

## Databases
The following databases are included in the *databases/* folder. Links go to the appropriate directory.
* [DroID](databases/DroID)
* [SignaLink2.0](databases/SignaLink)

## Statistics 

## Notes
- Ported over from bitbucket account at research/2016-06-17-apical-constriction/interactomes