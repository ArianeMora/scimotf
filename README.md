# scitf
[![codecov.io](https://codecov.io/github/ArianeMora/scitf/coverage.svg?branch=master)](https://codecov.io/github/ArianeMora/scimo?branch=master)
[![PyPI](https://img.shields.io/pypi/v/scitf)](https://pypi.org/project/scitf/)


sci-TF is a very simple package to help with finding motifs that are enriched in different clusters, that are also 
expressed in your dataset and make it easier to draw inferences on which TFs may be driving the observed changes.

The input to scitf is: 1) the output of fimo, fimo.tsv, 2) a csv file with gene identifier (e.g. name), cluster, log2FC,
 and p-value.

### Example format for fimo.tsv
``` 
motif_id        motif_alt_id    sequence_name   start   stop    strand  score   p-value q-value matched_sequence
SP5_MOUSE.H11MO.0.C             Gh      1668    1691    -       32.7879 9.78e-16        4e-09   GGGGGGGAGGGGGAGGGGGAGGGG
```

### Example format for cluster.csv
``` 
gene_name,cluster,log2FC,padj
Hoxa9,hindbrain,-2.8,0.00031
```

sci-TF will output two files, 1) scitf_detailed.csv, and 2) scitf_summary.csv. 

### scitf_detailed.csv

This gives a detailed output of each TF that was potentially able to bind to genes in a cluster.
``` 
cluster,motif,p-value,q-value,odds-ratio,count-genes-in-cluster,count-genes-bg,remainder-cluster,remainder-bg,tf-log2FC,tf-padj,tf-cluster,%-coverage,genes
```

### scitf_summary.csv
Gives a summary at the cluster level. i.e. what was the number of TFs that potentially bound to this cluster, and if 
any other clusters may be driving this one, each cluster will have 2 rows, one for TFs UP and one for TFs down.
``` 
cluster,motifs,unique-motifs,%-genes-with-motifs,%-genes-with-unique-motifs,direction,[list of other clusters with any TFs binding to found motifs]
```

### Overview
1) Filter fimo.tsv and remove any motifs that don't meet the p or qvalue threshold
2) Filter any motifs in fimo.tsv that don't exist in the users input data (have a 0 logFC)
3) For each TF for each cluster, count how many genes exist and perform a FET w.r.t the background
4) adjust p-values
5) summarise the identified TFs

Please post questions and issues related to sci-TF on the `Issues <https://github.com/ArianeMora/scitf/issues>`_  
section of the GitHub repository.
