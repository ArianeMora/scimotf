###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################
import pandas as pd
import numpy as np
import os
import sys
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from sciutil import SciUtil, SciException
from sciviso import Emapplot
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from wordcloud import WordCloud
import seaborn as sns
from tqdm import tqdm


class Epi2GeneException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class SciMotf_Doro:

    def __init__(self, doro_file: str, cluster_file: str, cluster_gene_id: str,
                 padj_protein: str, logfc_protein: str, padj_rna: str, logfc_rna: str,
                 output_dir='.', cluster_id='Regulation_Grouping_2', bg_cluster=None, cluster_pcutoff=1.0,
                 min_genes_in_cluster=3, tf_in_dataset=True, alpha=0.1, correction_method='fdr_bh',
                 min_odds_ratio=1.0, consider_direction=True):
        self.doro_file = doro_file
        self.consider_direction = consider_direction
        self.cluster_file = cluster_file
        self.tf_df, self.cluster_df = None, None
        self.c_id = cluster_id
        self.c_gid = cluster_gene_id
        self.padj_protein = padj_protein
        self.logfc_protein = logfc_protein
        self.padj_rna = padj_rna
        self.logfc_rna = logfc_rna
        self.logfc_protein = logfc_protein
        self.c_pcutoff = cluster_pcutoff
        self.output_dir = output_dir
        self.bg_cluster = bg_cluster
        self.bg_genes = []
        self.min_genes_in_cluster = min_genes_in_cluster
        self.tf_in_dataset = tf_in_dataset
        self.alpha = alpha
        self.correction_method = correction_method
        self.min_odds_ratio = min_odds_ratio
        self.u = SciUtil()

    def __load(self, doro_level):
        """ Check that the files exist and parse. """
        # Filter to doro level
        if doro_level is not None:
            self.tf_df = self.tf_df[self.tf_df['confidence'].isin(doro_level)]
        print(self.tf_df.head())

        # Before we run check that they aren't doing something dumb and that their output directory exists
        if not os.path.isdir(self.output_dir):
            raise Epi2GeneException(self.u.err_p(['Your output directory was not a directory., ', self.output_dir]))

        # Make a dataframe of the TFs and the genes they target
        reg_tfs = set(self.tf_df['tf'].values)
        tfs = []
        tf_targets = []
        tf_values = []
        target_values = []
        tf_padj = []
        protein_logFC_map = {}
        protein_sig_map = {}
        rna_logFC_map = {}
        rna_sig_map = {}
        # Check changes in transcription
        rna_logfc_values = self.cluster_df[self.logfc_rna].values
        rna_sig_values = self.cluster_df[self.padj_rna].values
        for i, gene_name in enumerate(self.cluster_df[self.c_gid].values):
            rna_logFC_map[gene_name] = rna_logfc_values[i]
            rna_sig_map[gene_name] = rna_sig_values[i]

        protein_logfc_values = self.cluster_df[self.logfc_protein].values
        protein_sig_values = self.cluster_df[self.padj_protein].values

        for i, gene_name in enumerate(self.cluster_df[self.c_gid].values):
            protein_logFC_map[gene_name] = protein_logfc_values[i]
            protein_sig_map[gene_name] = protein_sig_values[i]

        for tf in reg_tfs:
            targets = set(self.tf_df[self.tf_df['tf'] == tf]['target'].values)
            # FIrst check if we actually have a change by this TF
            if protein_sig_map.get(tf) and protein_sig_map.get(tf) < 1.0:
                for t in targets:
                    if rna_logFC_map.get(t) is not None and rna_sig_map.get(t) < 1.0:
                        # CHeck if the target is significantly changed IN the right direction
                        mor = self.tf_df[self.tf_df['tf'] == tf]
                        mor = mor[mor['target'] == t]['mor'].values[0]
                        if not self.consider_direction:
                            tfs.append(tf)
                            tf_targets.append(t)
                            tf_values.append(protein_logFC_map.get(tf))
                            target_values.append(protein_logFC_map.get(t))
                            tf_padj.append(protein_sig_map.get(tf))
                        else:
                            if mor == 1 and np.sign(rna_logFC_map.get(t)) == np.sign(protein_logFC_map.get(tf)):
                                tfs.append(tf)
                                tf_targets.append(t)
                                tf_values.append(protein_logFC_map.get(tf))
                                target_values.append(protein_logFC_map.get(t))
                                tf_padj.append(protein_sig_map.get(tf))
                            elif mor == -1 and np.sign(rna_logFC_map.get(t)) != np.sign(protein_logFC_map.get(tf)):
                                tfs.append(tf)
                                tf_targets.append(t)
                                tf_values.append(protein_logFC_map.get(tf))
                                target_values.append(protein_logFC_map.get(t))
                                tf_padj.append(protein_sig_map.get(tf))

        tf_val_df = pd.DataFrame()
        tf_val_df['motif_id'] = tfs
        tf_val_df['sequence_name'] = tf_targets
        tf_val_df['protein_TF_logFC'] = tf_values
        tf_val_df['rna_target_logFC'] = target_values
        tf_val_df['protein_TF_padj'] = tf_padj
        self.tf_df = tf_val_df

    def __gen_bb(self, rcm_clusters):
        self.cluster_df = pd.read_csv(self.cluster_file)
        self.tf_df = pd.read_csv(self.doro_file)
        """ Generate the background, if they don't pass a gene list with a background, then just use all the genes. """
        if self.bg_cluster:
            self.bg_genes = self.cluster_df[self.cluster_df[self.c_id] == self.bg_cluster][self.c_gid].values
            # Check the length is non-zero
            if len(self.bg_genes) == 0:
                msg = '\t'.join(['WARN: run: in generate background, you have no background genes?\n'
                                 ' Please check your column ID for the cluster:', str(self.c_id),
                                 '\n Also check your BG ID:', str(self.bg_cluster)])
                self.u.warn_p([msg])
                raise Epi2GeneException(msg)
        else:
            # Filter to include only the clusters we're intersted in
            tf_targets = list(set(self.tf_df['target'].values))
            if rcm_clusters is not None:
                self.bg_genes = list(set(tf_targets) & set(
                    self.cluster_df[self.cluster_df[self.c_id].isin(rcm_clusters)][self.c_gid].values))
                self.cluster_df = self.cluster_df[self.cluster_df[self.c_id].isin(rcm_clusters)]
            else:
                # Also need to add any targets in teh tf dataset that weren't included
                self.bg_genes = list(set(tf_targets) & set(self.cluster_df[self.c_gid].values))
            self.tf_df = self.tf_df[self.tf_df.target.isin(self.bg_genes)]
            self.cluster_df = self.cluster_df[self.cluster_df[self.c_gid].isin(self.bg_genes)]

            self.u.warn_p(['WARN: no background ID set, using all genes in the supplied DF as the background.\n'
                           'Number of genes: ', len(self.bg_genes)])

    def run(self, doro_level=['A'],
            rcm_clusters=["MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TPDE", "TPDE_TMDS", "TPDS", "TPDS_TMDE", ]):

        self.__gen_bb(rcm_clusters)
        self.__load(doro_level)
        # Keep track of the counts in each and also the genes that overlapped
        columns = ['Regulatory Cluster label', 'TF', 'p-value', 'odds-ratio',
                   'genes targeted by TF and in cluster',
                   'genes in cluster but not targeted by TF',
                   'genes in TF but not cluster',
                   'genes not in TF or cluster',
                   'TF padj',
                   'TF logFC',
                   'Mean logFC for genes targeted by TF and in cluster',
                   'Mean logFC for genes not targeted by TF and in cluster',
                   'gene_names']
        rows = []
        tf_cell_gene_ids = defaultdict(dict)
        gene_id = 'sequence_name'
        for regulatory_cluster in tqdm(self.cluster_df[self.c_id].unique()):
            r_df = self.cluster_df[self.cluster_df[self.c_id] == regulatory_cluster]
            # For each cell type, construct contingency table
            cluster_genes = r_df['external_gene_name'].values
            for tf in self.tf_df['motif_id'].unique():
                tf_df = self.tf_df[self.tf_df.motif_id == tf]
                tf_enriched_genes = tf_df[gene_id].values
                cont_table = np.zeros((2, 2))  # Rows=ct, not ct; Cols=Module, not Module
                if len(tf_enriched_genes) > 0:  # Some DE genes
                    in_tf_and_cluster = len(set(cluster_genes) & set(tf_enriched_genes))
                    cluster_not_tf = len(set(cluster_genes)) - in_tf_and_cluster
                    tf_not_cluster = len(set(tf_enriched_genes)) - in_tf_and_cluster
                    if in_tf_and_cluster > 2:  # Require there to be at least 3 genes overlapping
                        not_tf_not_cluster = len(set(self.bg_genes)) - (
                                    in_tf_and_cluster + cluster_not_tf + tf_not_cluster)
                        # Populating cont table
                        cont_table[0, 0] = in_tf_and_cluster
                        cont_table[1, 0] = cluster_not_tf
                        cont_table[0, 1] = tf_not_cluster
                        cont_table[1, 1] = not_tf_not_cluster
                        # Doing FET, Enrichment IN cell type, only.
                        odds_ratio, pval = fisher_exact(cont_table, alternative="greater")
                        genes_ids = list(set(cluster_genes) & set(tf_enriched_genes))
                        tf_cell_gene_ids[regulatory_cluster][tf] = genes_ids
                        # Lastly calculate the average RNA change for the genes targetted by the TF and
                        # also not tagrtted by the TF
                        tf_targted = list(set(cluster_genes) & set(tf_enriched_genes))
                        not_targeted = [c for c in cluster_genes if c not in tf_enriched_genes]
                        target_mean = np.mean(
                            self.cluster_df[self.cluster_df[self.c_gid].isin(tf_targted)][self.logfc_rna].values)
                        non_target_mean = np.mean(
                            self.cluster_df[self.cluster_df[self.c_gid].isin(not_targeted)][self.logfc_rna].values)

                        rows.append([regulatory_cluster, tf, pval, odds_ratio, in_tf_and_cluster,
                                     cluster_not_tf, tf_not_cluster, not_tf_not_cluster,
                                     tf_df['protein_TF_padj'].values[0],
                                     tf_df['protein_TF_logFC'].values[0],
                                     target_mean,
                                     non_target_mean,
                                     ' '.join(genes_ids)])

        odds_ratio_df = pd.DataFrame(data=rows, columns=columns)
        try:
            reg, padj, a, b = multipletests(odds_ratio_df['p-value'].values,
                                            alpha=0.05, method='fdr_bh', returnsorted=False)
            odds_ratio_df['p.adj'] = padj
        except:
            self.u.warn_p(["Unable to perform multiple tests, please run this yourself to adjust pvalues.",
                           "Returning unadjusted values... "])

        return odds_ratio_df


def plot_cluster_tf(filename, gene_ratio_min=1, padj_max=0.05, title='', fig_dir='',
                    label_font_size=9, figsize=(3, 3), axis_font_size=6,
                    rcm_labels=["MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TMDE", "TMDS", "TPDE", "TPDE_TMDS", "TPDS",
                                "TPDS_TMDE"],
                    save_fig=True):
    """

    Parameters
    ----------
    filename
    gene_ratio
    count_column
    padj
    overlap_column
    id_column
    label_column
    gene_ratio_min
    padj_max
    title
    label_font_size
    figsize
    axis_font_size
    min_count
    max_count
    min_overlap
    save_fig

    Returns
    -------

    """
    odds_ratio_df = pd.read_csv(filename)
    for r in rcm_labels:
        r_df = odds_ratio_df[odds_ratio_df['Regulatory Cluster label'] == r]
        r_df = r_df[r_df['genes targeted by TF and in cluster'] > gene_ratio_min]
        r_df = r_df[r_df['p.adj'] < padj_max]
        if len(r_df) > 1:
            eplot = Emapplot(r_df,
                             size_column='genes targeted by TF and in cluster',
                             color_column='p.adj',
                             id_column='TF',
                             label_column='TF',
                             overlap_column='gene_names', overlap_sep=' ', title=r,
                             config={'figsize': figsize, 'label_font_size': label_font_size,
                                     'axis_font_size': axis_font_size})
            eplot.build_graph()
            plt.title(title)
            plt.gca().set_clip_on = False
            if save_fig:
                plt.savefig(f'{fig_dir}TF_{title.replace(" ", "-")}_network.svg', bbox_inches='tight',
                            transparent=True)
            plt.show()

            x, y = np.ogrid[:300, :300]

            mask = (x - 150) ** 2 + (y - 150) ** 2 > 130 ** 2
            mask = 255 * mask.astype(int)
            wordfeqs = defaultdict(int)
            for g in r_df['gene_names'].values:
                for w in g.split(' '):
                    w = w.replace(' ', '.')
                    wordfeqs[w] += 1
            total_words = len(wordfeqs)
            for w in wordfeqs:
                wordfeqs[w] = wordfeqs[w] / total_words
            wordcloud = WordCloud(background_color="white", mask=mask, colormap='viridis',
                                  repeat=False).generate_from_frequencies(wordfeqs)

            plt.figure()
            plt.rcParams['svg.fonttype'] = 'none'  # Ensure text is saved as text
            plt.rcParams['figure.figsize'] = figsize
            font_family = 'sans-serif'
            font = 'Arial'
            sns.set(rc={'figure.figsize': figsize, 'font.family': font_family,
                        'font.sans-serif': font, 'font.size': 12}, style='ticks')
            plt.figure()
            plt.imshow(wordcloud, interpolation="bilinear")
            plt.axis("off")
            if save_fig:
                wordcloud_svg = wordcloud.to_svg(embed_font=True)
                f = open(f'{fig_dir}{title}_{r}_WordCloud.svg', "w+")
                f.write(wordcloud_svg)
                f.close()
                plt.savefig(f'{fig_dir}{title}_{r}_WordCloud.png', bbox_inches='tight')
            plt.show()