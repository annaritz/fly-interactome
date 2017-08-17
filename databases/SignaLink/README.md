# [SignaLink2.0: A signaling pathway resource with multi-layered regulatory networks(http://signalink.org)
- [SignaLink 2.0 - A signaling pathway resource with multi-layered regulatory networks.](http://www.biomedcentral.com/1752-0509/7/7) Fazekas D, Koltai M, Türei D, Módos D, Pálfy M, Dúl Z, Zsákai L, Szalay-Bekő M, Lenti K, Farkas I J, Vellai T, Csermely P, Korcsmáros T. BMC Systems Biology 2013, 7:7.
- [Uniformly curated signaling pathways reveal tissue-specific cross-talks and support drug target discovery.](http://bioinformatics.oxfordjournals.org/cgi/reprint/btq310?ijkey=exmTiN0PbYl8mPH&keytype=ref) Korcsmáros T, Farkas I J, Szalay M S, Rovó P, Fazekas D, Spiro Z, Böde C, Lenti K, Vellai T, Csermely P. Bioinformatics 26:2042-2050 (2010).

Downloaded 2016-07-27 from the [Download Page](http://signalink.org/download) with the following requirements:
- **Species:** Drosophila
- **Interaction Types:**  Further Interactions, Directed protein-protein interactions, Post-translational modficators, Pathway regulators, Pathway Members (everything except Transcriptional Regulators)
-**Pathways:** All
-**TF Network:** Do not include
-**Output Format:** csv
The output file is `signalink-download.csv`

FlyBase-mapped interactome is `signalink-flybase.txt` and UniProt-mapped interactome is `signalink-uniprot.txt`.

## Running the SignaLink Parser

```python build-signalink.py```

The file contains both FlyBase IDs and UniProt IDs.  The code builds interactomes with 5236 uniprot edges and 2343 flybase edges, respectively.