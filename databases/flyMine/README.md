# [FlyMine: An integrated database for Drosophila genomics](http://www.flymine.org/)
- Lyne R, Smith R, Rutherford K, Wakeling M, Varley A, Guillier F, Janssens H, Ji W, Mclaren P, North P, Rana D, Riley T, Sullivan J, Watkins X, Woodbridge M, Lilley K, Russell S, Ashburner M, Mizuguchi K, Micklem G. [FlyMine: an integrated database for Drosophila and Anopheles genomics.](https://www.ncbi.nlm.nih.gov/pubmed/17615057) Genome Biol. 2007;8(7):R129.

FlyBase-mapped interactome is `flymine-flybase.txt`.  UniProt-mapped interactome is `flymine-uniprot.txt`.

## Get the data


Built a query to get all interactions.  Used the [QueryBuilder](http://www.flymine.org/flymine/customQuery.do). Select "Interaction" as the data type to begin the query, and constraint to be Interaction.Type = physical.  Select the following fields to show in the results:
- Participant 1 (Name, DBidentifier, Secondary identifier, Symbol)
- Participant 2 (Name, DBidentifier, Secondary identifier, Symbol)
- Interaction (Name, Experiment.publication.pubmedID, Experiment.InteractionDetectionMethods.Name)
Select the Python link and copy/paste the code into a file.  This is stored here as `auto-gen-script.py,` with some extra code at the bottom to write the output to `flymine-query.txt`.  It requires [intermine](http://intermine.readthedocs.org/en/latest/web-services/).

```python auto-gen-script.py```

## Running the flyMine parser

Once you have the interactions in `flymine-query.txt`, you can run the parser:

```
python build-flymine.py
```

Parses 278370 edges.

## Convert the FlyBase Interactome into the UniProt Interactome

```
python ../../utils/map.py --mapfrom FLYBASE --mapto uniprot flymine-flybase.txt flymine-uniprot.txt 
```
The output:
```
278370 lines in original file flymine-flybase.txt
225154 lines written to flymine-uniprot.txt
180045 lines were missing an id
```

