# Fly Interactome Script and Data Files

The files include the following:

- `interactome-uniprot.txt`: Interactome in [UniProtKB IDs](http://www.uniprot.org/)
- `interactomeflybase.txt`: Interactome in [FlyBase IDs](http://flybase.org/).

The namespace mapping occurs at the specific database level - see the scripts in the [databases/](../databases/) directory.  

## Building the Interactome

```
python build-fly-interactome.py
```

- TODO: add argument parsing for this script.
- TODO: Right now it assumes that all interactions are undirected. This is a simplistic assumption.
- TODO: add summary statistics.