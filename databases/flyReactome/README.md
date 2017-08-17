# [FlyReactome: a curated knowledgebase of drosophila melanogaster pathways](http://fly.reactome.org/)
- Matthews L, Gopinath G, Gillespie M, Caudy M, Croft D, de Bono B, Garapati P, Hemish J, Hermjakob H, Jassal B, Kanapin A, Lewis S, Mahajan S, May B, Schmidt E, Vastrik I, Wu G, Birney E, Stein L, D'Eustachio P. [Reactome knowledgebase of biological pathways and processes.](http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=18981052) Nucleic Acids Res. 2008 Nov 3.

Downloaded 2016-07-27 from the [Download Page](http://fly.reactome.org/download/index.html). Selected the file "Fly protein-protein interaction pairs in tab-delimited format." Take note of the [inferred reaction details](http://fly.reactome.org/download/interactions.README.txt) from the page.  FlyBase-mapped interactome is `flyReactome-flybase.txt`.  UniProt-mapped interactome is `flyReactome-uniprot.txt`.

## Running the FlyReactome Parser

Download the file from the [Download Page](http://fly.reactome.org/download/index.html), unzip it and place it in the same directory as `build-flyReactome.py`.  The file should be named `drosophila_melanogaster.interactions.txt`.

```
python build-flyReactome.py
```

Writes 612 interactions to the uniprot file.

## Convert the UniProt interactome to the FlyBase interactome

```
python ../../utils/map.py --mapfrom uniprot --mapto FLYBASE flyReactome-uniprot.txt flyReactome-flybase.txt
```

The output:
```
612 lines in original file flyReactome-uniprot.txt
468 lines written to flyReactome-flybase.txt
144 lines were missing an id
```