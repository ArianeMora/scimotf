# sci-moTF
[![PyPI](https://img.shields.io/pypi/v/scimotf)](https://pypi.org/project/scitf/)

sci-moTF is a very simple package to help with finding motifs that are enriched in different clusters, that are also 
expressed in your dataset and make it easier to draw inferences on which TFs may be driving the observed changes.

## Install

```
pip install scimotf
```

There are two ways to run scimotf, 1) using DoRoTHea, 2) using FIMO.

## Example using [DoRothEA](https://bioconductor.org/packages/release/data/experiment/html/dorothea.html):

```
from scimotf import SciMotf_Doro

rcm_file = f'file output from the scircm package'
tf_file = 'dorothea_hs_ABCD.csv' # File downloaded from DoRothEA
mo = SciMotf_Doro(doro_file=tf_file, cluster_file=rcm_file, 
                 cluster_gene_id='external_gene_name', # got to match motif
                 padj_protein='column with your protein padj value',
                  logfc_protein='column with the protein logFC', 
                  padj_rna='column with the RNA padj',
                  logfc_rna='column with the RNA logFC', 
                  output_dir='')

# Run with the letters your interested in (i.e. A, B, C, D) see doro paper for deets
df = mo.run(['A'], rcm_clusters=["TMDE", "TMDS", "MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TPDE", "TPDE_TMDS", "TPDS", "TPDS_TMDE",])
df.to_csv(f'scimotif_DORO_A.csv')
```

#### Plot the results

```
from scimotf import plot_cluster_tf
plot_cluster_tf(f'scimotif_DORO_A.csv', save_fig=True, fig_dir='')
```

## Example using FIMO:
The input to scimotf is: 1) the output of [FIMO](https://meme-suite.org/meme/doc/fimo.html?man_type=web>) , fimo.tsv, 2) a csv file with gene identifier (e.g. name), cluster, log2FC,
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

### sictf_motif_merged_fp-0.05_cp-1.0.csv

This gives a detailed output of each TF that was potentially able to bind to genes in a cluster.
``` 
cluster,motif,p-value,q-value,odds-ratio,count-genes-in-cluster,count-genes-bg,remainder-cluster,remainder-bg,tf-log2FC,tf-padj,tf-cluster,%-coverage,genes
```

### Overview
1) Filter fimo.tsv and remove any motifs that don't meet the p or qvalue threshold
2) Filter any motifs in fimo.tsv that don't exist in the users input data (have a 0 logFC)
3) For each TF for each cluster, count how many genes exist and perform a FET w.r.t the background
4) adjust p-values
5) summarise the identified TFs

Please post questions and issues related to sci-moTF on the [Issues page](https://github.com/ArianeMora/scimotf/issues)_  

section of the GitHub repository.

