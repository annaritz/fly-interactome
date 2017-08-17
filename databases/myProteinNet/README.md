# [MyProteinNet: Build Up-to-date protein interaction networks for organisms, tissues and user-defined contexts](http://netbio.bgu.ac.il/myproteinnet/)
- Omer Basha, Dvir Flom, Ruth Barshir, Ilan Smoly, Shoval Tirman & Esti Yeger-Lotem.
[MyProteinNet: Build up-to-date protein interaction networks for organisms, tissues and user-defined contexts.](http://nar.oxfordjournals.org/content/43/W1/W258.full) Nucl. Acids Res. (1 July 2015) 43 (W1): W258-W263.

You can submit a job on the [main web page](http://netbio.bgu.ac.il/myproteinnet/index.html).  There was a bug in biomart when I downloaded this version in 2016-08-15 and Omar Basha generated the interactome.  The files here are copied directly from that job.
- `GlobalInteractome.tsv`: The unfiltered interactome in ensmbel identifiers (I think these are FlyBase identifiers, actually)
- `GlobalInteractomeGeneSymbol.tsv`: The unfiltered interactome in gene symbol identifiers.
- `InteractionsDetectionMethods.tsv`: This file contains the interactions with their detection methods in the MITAB format.

 FlyBase-mapped interactome is `myProteinNet-flybase.txt`.  UniProt-mapped interactome is `myProteinNet-uniprot.txt`.  There are no PubMedIDs included in the output files (all are listed as "NA"). 

 ## Running the myProteinNet parser

 The parser relies on converting MI identifiers to common names. Thus, I downloaded the Molecular Interaction (MI) terms from the [Ontology Lookup Service](https://www.ebi.ac.uk/ols/ontologies/mi).  This file is stored as `mi.owl`, and I rely on the fact that the MI term is followed immediately by the name.

```
python build-myproteinnet.py
```

1489 MI terms mapped and 41530 edges read.

## Convert the FlyBase Interactome into the UniProt Interactome

```
python ../../utils/map.py --mapfrom FLYBASE --mapto uniprot myProteinNet-flybase.txt myProteinNet-uniprot.txt 
```

The output:

```
41530 lines in original file myProteinNet-flybase.txt
87529 lines written to myProteinNet-uniprot.txt
477 lines were missing an id
```
