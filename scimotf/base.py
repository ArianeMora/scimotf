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
from scipy import stats
from statsmodels.stats.multitest import multipletests

from sciutil import SciUtil, SciException


class Epi2GeneException(SciException):
    def __init__(self, message=''):
        Exception.__init__(self, message)


class SciMotf:

    def __init__(self, fimo_file: str, cluster_file: str, cluster_id: str, cluster_gene_id: str,
                 cluster_p: str, cluster_logfc: str, output_dir: str, bg_cluster=None, fimo_pcol='q-value',
                 fimo_pcutoff=0.05, cluster_pcutoff=1.0, min_genes_in_cluster=3,
                 tf_in_dataset=True, alpha=0.1, correction_method='fdr_bh', min_odds_ratio=2.0):
        self.fimo_file = fimo_file
        self.cluster_file = cluster_file
        self.fimo_df, self.cluster_df = None, None
        self.c_id = cluster_id
        self.c_gid = cluster_gene_id
        self.c_p = cluster_p
        self.c_logfc = cluster_logfc
        self.f_p = fimo_pcol
        self.f_pcutoff = fimo_pcutoff
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

    def __load(self):
        """ Check that the files exist and parse. """
        self.fimo_df = pd.read_csv(self.fimo_file, sep='\t')
        print(self.fimo_df.head())
        self.cluster_df = pd.read_csv(self.cluster_file)
        # Before we run check that they aren't doing something dumb and that their output directory exists
        if not os.path.isdir(self.output_dir):
            raise Epi2GeneException(self.u.err_p(['Your output directory was not a directory., ', self.output_dir]))
        # Otherwise group on the motif ID and filter on pvalue (or qvalue)
        fimo = pd.DataFrame(self.fimo_df[self.fimo_df[self.f_p] < self.f_pcutoff])
        # Group by motifs and have a look at how many there are
        fimo['count'] = np.ones(len(fimo))
        fimo['g_id'] = fimo['motif_id'].values
        fimo_grped = fimo.set_index('motif_id')
        fimo_grped = fimo_grped.groupby('motif_id')
        self.u.dp(['Fimo raw: ', len(fimo), 'Fimo grouped: ', len(fimo_grped)])
        self.fimo_df = fimo_grped

    def __gen_bb(self):
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
            self.bg_genes = self.cluster_df[self.c_gid].values
            self.u.warn_p(['WARN: no background ID set, using all genes in the supplied DF as the background.\n'
                           'Number of genes: ', len(self.bg_genes)])

    def __calc_adj_matrix(self):
        num_bg = len(self.bg_genes)
        cluster_totals = {'bg': num_bg}
        # Get the unique cluster groups
        clusters = set(list(self.cluster_df[self.c_id].values))
        cluster_genes = {}
        for i, cluster in enumerate(clusters):
            if i != self.bg_cluster and i is not None and i != "None":  # Potentially empty group.
                # Get the number of genes in each cluster
                cluster_genes[cluster] = set(list(self.cluster_df[self.cluster_df[self.c_id] == cluster][self.c_gid].values))
                cluster_totals[f'{cluster}'] = len(cluster_genes[cluster])

        gene_dict = {'bg': {}}
        motif_dict_reg_grps = {'bg': {}}
        all_gene_dict = {}
        for g in self.bg_genes:
            gene_dict['bg'][g.lower()] = True
            all_gene_dict[g.lower()] = g

        for i, cluster in enumerate(cluster_genes):
            gene_dict[f'{cluster}'] = {}
            motif_dict_reg_grps[f'{cluster}'] = {}
            for g in cluster_genes[cluster]:
                gene_dict[f'{cluster}'][g.lower()] = True
                all_gene_dict[g.lower()] = g

        # Save the results to a file
        tf_info = {}
        motif_to_genes = {}

        for motif in self.fimo_df:
            m = motif[0].split('_')[0].lower().split('_')[0]
            motif_id = motif[0]
            # First check if this motif exists as a TF in our dataset if the user has selected this
            # ToDo: unit test for this
            if not self.tf_in_dataset or (self.tf_in_dataset and all_gene_dict.get(m)):
                # let's add the information for this tf (we want 1) log2FC, padj and cluster this belongs to
                tf = self.cluster_df[self.cluster_df[self.c_gid] == all_gene_dict.get(m)]
                if len(tf) > 0:  # If this was in our dataset
                    tf_info[motif_id] = {'log2fc': tf[self.c_logfc].values[0], 'padj': tf[self.c_p].values[0],
                                  'cluster': tf[self.c_id].values[0]}
                motif_genes = motif[1]['sequence_name'].values
                # split the sequence genes to be the final name ToDo Modularise
                motif_genes = [ms.split('|')[-1] for ms in motif_genes]
                # Keep track of the genes associated with this motif
                motif_to_genes[motif_id] = list(motif_genes)
                for i in motif_dict_reg_grps:
                    motif_dict_reg_grps[i][motif_id] = []
                for g in motif_genes:
                    for k, cluster_gene_map in gene_dict.items():
                        if cluster_gene_map.get(g.lower()):
                            motif_dict_reg_grps[k][motif_id].append(all_gene_dict.get(g.lower()))
                # Lastly let's turn it into a set
                motif_dict_reg_grps['bg'][motif_id] = set(motif_dict_reg_grps['bg'][motif_id])
                for k, v in gene_dict.items():
                    motif_dict_reg_grps[k][motif_id] = set(motif_dict_reg_grps[k][motif_id])
        return motif_dict_reg_grps, motif_to_genes, cluster_totals, tf_info

    def __gen_output(self, motif_dict_reg_grps: dict, motif_to_genes: dict, cluster_totals: dict, tf_info: dict):
        """ Make summary and TF file. """
        num_bg_all = cluster_totals['bg']
        fout_name = f'{self.output_dir}sictf_motif_merged_fp-{self.f_pcutoff}_cp-{self.c_pcutoff}.tsv'

        with open(fout_name, 'w+') as m_out:
            # cluster,motif,p-value,q-value,odds-ratio,count-genes-in-cluster,count-genes-bg,remainder-cluster,remainder-bg,
            # tf-log2FC,tf-padj,tf-cluster,%-coverage,genes

            m_out.write(
                'cluster\tmotif\tp-value\tq-value\todds-ratio\tcount-genes-in-cluster\tcount-genes-bg'
                '\tremainder-cluster\tremainder-bg\ttf-log2FC\ttf-padj\ttf-cluster\t%-coverage\tgenes\n')

            for m in motif_to_genes:

                bg_genes = motif_dict_reg_grps['bg'][m]
                num_bg_m = len(bg_genes)

                for k, v in motif_dict_reg_grps.items():

                    if k != 'bg' and k != "None" and k is not None:
                        if len(motif_dict_reg_grps[k][m]) >= self.min_genes_in_cluster:
                            num_in_cluster = len(motif_dict_reg_grps[k][m])
                            if not self.bg_cluster:
                                # BG also contains that cluster so we need to subtract this out
                                num_in_bg = num_bg_m - num_in_cluster
                                num_bg = num_bg_all - cluster_totals[k]
                            else:
                                # Otherwise it is assumed they are independent
                                num_bg = num_bg_all
                            oddsratio, pvalue = stats.fisher_exact([[num_in_cluster, num_in_bg],
                                                                    [cluster_totals[k] - num_in_cluster,
                                                                     num_bg]])
                            if pvalue < self.c_pcutoff:
                                if num_in_cluster == cluster_totals[k]:
                                    perc = 100
                                else:
                                    perc = 100 * (1.0 * num_in_cluster/(cluster_totals[k]))
                                m_tf_info = tf_info.get(m)
                                # Filter by a minimum enrichment
                                if oddsratio > self.min_odds_ratio:
                                    if m_tf_info:
                                        m_out.write(f'{k}\t{m}\t{pvalue}\t{0}\t{oddsratio}\t{num_in_cluster}\t{num_in_bg}\t'
                                                    f'{cluster_totals[k] - num_in_cluster}\t{num_bg - num_in_bg}\t'
                                                    f'{m_tf_info["log2fc"]}\t{m_tf_info["padj"]}\t{m_tf_info["cluster"]}\t'
                                                    f'{perc}\t'
                                                    f'{",".join(motif_dict_reg_grps[k][m])}\n')
                                    else:
                                        m_out.write(f'{k}\t{m}\t{pvalue}\t{0}\t{oddsratio}\t{num_in_cluster}\t{num_in_bg}\t'
                                                    f'{cluster_totals[k] - num_in_cluster}\t{num_bg - num_in_bg}\t'
                                                    f'{None}\t{None}\t{None}\t'
                                                    f'{perc}\t'
                                                    f'{",".join(motif_dict_reg_grps[k][m])}\n')

        # Re-read in the file to perform correction for FDR and to generate teh summary file
        # https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
        m_df = pd.read_csv(fout_name, sep='\t')
        reg, padj, a, b = multipletests(m_df['p-value'].values, alpha=self.alpha, method=self.correction_method,
                                        returnsorted=False)
        m_df['q-value'] = padj
        m_df.to_csv(fout_name, sep='\t', index=False)
        self.motif_df = m_df
        return m_df

    def add_tf_predictions(self, motif_df: pd.DataFrame, prediction_file: str,
                                 tf_column: str, value_column: str, column_lbl: str,
                                 output_filename=None):
        """ This can be done after __gen_output """
        if self.motif_df is None:
            self.u.err_p(["ERR: You need to first run scimo before you can add the predictions."])
            return
        # Otherwise we want to iterate through the file of TFs and basically just match the TF to the TF in the motif
        # file and add in the value column.
        # First make a dictionary mapping the tf to the value from dorethea
        tf_df = pd.read_csv(prediction_file)
        tf_to_value = {}
        values = tf_df[value_column].values
        for i, c in enumerate(tf_df[tf_column].values):
            tf_to_value[c.lower()] = values[i]

        # Now assign these to the motif dataframe
        motif_values_lst = []
        for m in motif_df['motif'].values:
            motif_values_lst.append(tf_to_value.get(m.lower().split('_')[0]))
        # add as a column
        motif_d_df = motif_df.copy()
        motif_d_df[f'{column_lbl}_{value_column}'] = motif_values_lst
        if output_filename:
            # Save
            motif_d_df.to_csv(output_filename, sep='\t', index=False)
        return motif_d_df

    def run(self):
        self.__load()
        self.__gen_bb()
        # Now we want to perform the FET on each of the clusters
        motif_dict_reg_grps, motif_to_genes, cluster_totals, tf_info = self.__calc_adj_matrix()
        m_df = self.__gen_output(motif_dict_reg_grps, motif_to_genes, cluster_totals, tf_info)
        return m_df
