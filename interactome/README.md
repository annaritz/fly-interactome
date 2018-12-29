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

## Map Evidence Codes to MI Ontology

Get the [Molecular Interactions Controlled Vocabulary](https://www.ebi.ac.uk/ols/ontologies/mi) (click the Download button).

```
python ../utils/evidence-to-mi.py mi.owl interactome-flybase.txt interactome-flybase-evidence.txt 4
python ../utils/evidence-to-mi.py mi.owl interactome-flybase-collapsed.txt interactome-flybase-collapsed-evidence.txt 6
python ../utils/evidence-to-mi.py mi.owl interactome-uniprot.txt interactome-uniprot-evidence.txt 4
python ../utils/evidence-to-mi.py mi.owl interactome-uniprot-collapsed.txt interactome-uniprot-collapsed-evidence.txt 6
```

For example, merging IDs in the Flybase collapsed interactome (`interactome-flybase-collapsed.txt`) maps these types of evidence codes:
```
81 evidence types
set(['MI:0416', 'worm_interologs', 'Post-translational_modification', 'MI:0091', 'MI:0402', 'human_interologs', 'indirect_complex', 'MI:0428', 'MI:0405', 'MI:0404', 'MI:0424', 'MI:0435', 'curagen_yth', 'MI:0423', 'MI:0813', 'MI:0065', 'MI:0067', 'MI:0892', 'MI:0040', 'MI:0588', 'MI:0084', 'Directed_protein-protein_interaction', 'MI:0676', 'MI:0007', 'MI:0006', 'MI:0049', 'MI:0004', 'MI:0809', 'yeast_interologs', 'MI:0027', 'MI:0029', 'MI:0028', 'interaction_from_external_databases', 'MI:0047', 'MI:0226', 'MI:0045', 'FB2015_05', 'dpim_coapcomplex', 'MI:0963', 'fly_other_physical', 'flybase_ppi', 'MI:0077', 'MI:0410', 'Interaction_between_pathway_members', 'MI:0413', 'MI:0727', 'perrimon_coapcomplex', 'MI:0107', 'direct_complex', 'MI:0729', 'hub-spoke-model', 'MI:0434', 'MI:0686', 'reaction', 'MI:0071', 'hybrigenics_yth', 'finley_yth', 'MI:0826', 'MI:0411', 'MI:0663', 'MI:0018', 'MI:0096', 'MI:0081', 'MI:0090', 'MI:0401', 'MI:0889', 'MI:0055', 'MI:0397', 'MI:0030', 'MI:0276', 'MI:0400', 'MI:0114', 'MI:0013', 'MI:0406', 'MI:0051', 'MI:0254', 'MI:0053', 'anti-tag-coip', 'MI:0019', 'MI:0399', 'MI:0415'])
58 mapped types
   filter-binding MI:0049
   three-hybrid MI:0588
   far-western-blotting MI:0047
   two-hybrid-array MI:0397
   bimolecular-fluorescence-complementation MI:0809
   two-hybrid MI:0018
   fluorescence-microscopy MI:0416
   phage-display MI:0084
   protein-complementation-assay MI:0090
   isothermal-titration-calorimetry MI:0065
   two-hybrid-fragment-pooling-approach MI:0399
   cosedimentation-through-density-gradient MI:0029
   enzymatic-study MI:0415
   acetylase-assay MI:0889
   genetic-interference MI:0254
   fluorescent-resonance-energy-transfer MI:0055
   comigration-in-non-denaturing-gel-electrophoresis MI:0404
   tandem-affinity-purification MI:0676
   surface-plasmon-resonance MI:0107
   confocal-microscopy MI:0663
   biophysical MI:0013
   x-ray-crystallography MI:0114
   affinity-chromatography-technology MI:0004
   cosedimentation MI:0027
   deacetylase-assay MI:0406
   fluorescence-technology MI:0051
   nuclear-magnetic-resonance MI:0077
   cross-linking-study MI:0030
   anti-bait-coimmunoprecipitation MI:0006
   protease-assay MI:0435
   fluorescence-polarization-spectroscopy MI:0053
   experimental-interaction-detection MI:0045
   electrophoretic-mobility-shift-assay MI:0413
   interactome-parallel-affinity-capture MI:0963
   molecular-sieving MI:0071
   blue-native-page MI:0276
   coimmunoprecipitation MI:0019
   electron-tomography MI:0410
   anti-tag-coimmunoprecipitation MI:0007
   chromatin-immunoprecipitation-assay MI:0402
   cosedimentation-in-solution MI:0028
   solid-phase-assay MI:0892
   unspecified-method MI:0686
   enzyme-linked-immunosorbent-assay MI:0411
   proximity-ligation-assay MI:0813
   luminescence-based-mammalian-interactome-mapping MI:0729
   pull-down MI:0096
   ion-exchange-chromatography MI:0226
   chromatography-technology MI:0091
   protein-kinase-assay MI:0424
   affinity-technology MI:0400
   phosphatase-assay MI:0434
   in-gel-kinase-assay MI:0423
   light-scattering MI:0067
   imaging-technique MI:0428
   electron-microscopy MI:0040
   biochemical MI:0401
   peptide-array MI:0081

```

## Developer Notes

- TODO: add argument parsing for this script.
- TODO: Right now it assumes that all interactions are undirected. This is a simplistic assumption.
- TODO: add summary statistics.
- TODO: add common names.