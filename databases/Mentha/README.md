# [Mentha: The Interactome Browser](http://mentha.uniroma2.it)
- [Mentha: a resource for browsing integrated protein-interaction networks.](http://dx.doi.org/10.1038/nmeth.2561)
Alberto Calderone, Luisa Castagnoli & Gianni Cesareni
Nature Methods 10, 690 (2013).

Downloaded 2017-08-16 from the [Download Page](http://mentha.uniroma2.it/download.php); selected *Drosophila melanogaster* only.  FlyBase-mapped interactome is `mentha-flybase.txt`.  UniProt-mapped interactome is `mentha-uniprot.txt`.

## Running the Mentha Parser

Download the Mentha *Drosophila* interactome and unzip it.  The file should be named `7227` and reside in the same directory as `build-mentha.py`.  It only contains UniProt IDs.

```python build-mentha.py```

Writes 45669 lines to mentha-uniprot.txt. 

## Convert Uniprot Interactome to the FlyBase Interactome

```
python ../../utils/map.py --mapfrom uniprot --mapto FLYBASE mentha-uniprot.txt mentha-flybase.txt 
```

The output:
```
45669 lines in original file mentha-uniprot.txt
43514 lines written to mentha-flybase.txt
3525 lines were missing an id
```