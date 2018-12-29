# [DroID: The Comprehensive Drosophila Interactions Database](http://www.droidb.org/)
* Murali T, Pacifico S, Yu J, Guest S, Roberts GG 3rd, Finley RL Jr. [DroID 2011: a comprehensive, integrated resource for protein, transcription factor, RNA and gene interactions for Drosophila.](http://nar.oxfordjournals.org/content/early/2010/10/28/nar.gkq1092.long) Nucleic Acids Res. 2010 Oct 29.
* Yu J, Pacifico S, Liu G, Finley RL Jr. [DroID: the Drosophila Interactions Database, a comprehensive resource for annotated gene and protein interactions.](http://www.biomedcentral.com/1471-2164/9/461) BMC Genomics. 2008 Oct 7;9:461.

Downloaded 2016-07-25 from the [Download Page](http://www.droidb.org/Downloads.jsp) (version 2015_12). FlyBase-mapped interactome is `DroID-flybase.txt`.  UniProt-mapped interactome is `DroID-uniprot.txt`.

## Running the DroID parser

Download DroID files from the [Download Page](http://www.droidb.org/Downloads.jsp) and unzip it in the directory where `build-DroID.py` resides.

```
python build-DroID.py 
```

The output:
```
curagen_yth has 19507 lines
finley_yth has 9008 lines
hybrigenics_yth has 1843 lines
fly_other_physical has 18238 lines
flybase_ppi has 19478 lines
dpim_coapcomplex has 61260 lines
perrimon_coapcomplex has 385 lines
human_interologs has 68786 lines
worm_interologs has 3650 lines
yeast_interologs has 81791 lines
Wrote 262179 lines to DroID.txt (some edges may be supported by different evidence sources)

316 were missing PubMed IDs and were replaced by NAs.
```

## Convert the FlyBase Interactome into the UniProt Interactome

```
python ../../utils/map.py --mapfrom FLYBASE --mapto uniprot DroID-flybase.txt DroID-uniprot.txt 
```
The output (as of Dec 2018):
```
262437 lines in original file DroID-flybase.txt
496848 lines written to DroID-uniprot.txt
6111 lines were missing an id
```
