# Fly Interactome Script and Data Files

The files include the following:

- `interactome-uniprot.txt`: Interactome in [UniProtKB IDs](http://www.uniprot.org/)
- `interactomeflybase.txt`: Interactome in [FlyBase IDs](http://flybase.org/).

The namespace mapping occurs at the specific database level - see the scripts in the [databases/](../databases/) directory.  Files containing the mapped namespaces for each interactome are `nodes-uniprot.txt` and `nodes-flybase.txt`, whicy use the `mygene` Python module for mapping (see the scripts in the [../utils/](../utils/) directory).

## Building the Interactome

```
python build-fly-interactome.py
```

## Get Mapping of Nodes


```
python ../utils/map-interactome-nodes.py --mapfrom uniprot --mapto FLYBASE interactome-uniprot.txt nodes-uniprot.txt
python ../utils/map-interactome-nodes.py --mapfrom FLYBASE --mapto uniprot interactome-flybase.txt nodes-flybase.txt
```

## Collapse Interactome

Map the interactome to common name (merging names if they correspond to the same identifier).  Two edges that have the same nodes are collapsed into one line, with evidence sources concatenated.

```
python ../utils/map.py --mapfrom uniprot --mapto symbol interactome-uniprot.txt interactome-uniprot-collapsed.txt --merge_ids
python ../utils/map.py --mapfrom FLYBASE --mapto symbol interactome-flybase.txt interactome-flybase-collapsed.txt --merge_ids
```

## Developer Notes

- TODO: add argument parsing for this script.
- TODO: Right now it assumes that all interactions are undirected. This is a simplistic assumption.
- TODO: add summary statistics.
- TODO: add common names.