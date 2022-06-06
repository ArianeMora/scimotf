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

import os
import shutil
import tempfile
import unittest

from scimotf import SciMotf


class TestClass(unittest.TestCase):

    @classmethod
    def setup_class(self):
        local = True
        # Create a base object since it will be the same for all the tests
        THIS_DIR = os.path.dirname(os.path.abspath(__file__))

        self.data_dir = os.path.join(THIS_DIR, 'data/')
        if local:
            self.tmp_dir = os.path.join(THIS_DIR, 'data/tmp/')
            if os.path.exists(self.tmp_dir):
                shutil.rmtree(self.tmp_dir)
            os.mkdir(self.tmp_dir)
        else:
            self.tmp_dir = tempfile.mkdtemp(prefix='scimo_tmp_')
        # Setup the default data for each of the tests
        self.fimo = os.path.join(self.data_dir, 'fimo_sml.tsv')
        self.fimo_all = os.path.join(self.data_dir, 'fimo.tsv')
        self.cluster = os.path.join(self.data_dir, 'cluster.csv')
        self.mm10_annot = os.path.join(self.data_dir, 'mmusculus_gene_ensembl-GRCm38.p6.csv')
        self.hg38_annot = os.path.join(self.data_dir, 'hsapiens_gene_ensembl-GRCh38.p13.csv')

    @classmethod
    def teardown_class(self):
        shutil.rmtree(self.tmp_dir)


class TestScimotf(TestClass):

    def test_scimotf(self):
        self.setup_class()
        """ Tests the generic function of scimo """
        # fimo_file: str, cluster_file: str, cluster_id: str, cluster_gene_id: str,
        #                  cluster_p: str, cluster_logfc: str, output_dir: str
        scimo = SciMotf(self.fimo, self.cluster, "cluster", "gene_name", "padj", "log2FC", self.data_dir,
                        tf_in_dataset=False)
        m_df = scimo.run()
        # Now we want to assert some simple things since we have validated the output
        assert m_df[m_df['motif'] == 'MAZ_MOUSE.H11MO.0.C']['count-genes-in-cluster'].values[0] == 2
        assert m_df[m_df['motif'] == 'WT1_MOUSE.H11MO.0.C']['remainder-bg'].values[0] == 10
        assert len(m_df[m_df['cluster'] == 'a']) == 3

    # def test_scimotf_SIRCLE(self):
    #     self.setup_class()
    #     # fimo_file: str, cluster_file: str, cluster_id: str, cluster_gene_id: str,
    #     # cluster_p: str, cluster_logfc: str, output_dir: str, bg_cluster = None, fimo_pcol = 'q-value',
    #     # fimo_pcutoff = 0.05, cluster_pcutoff = 1.0, min_genes_in_cluster = 3,
    #     # tf_in_dataset = True, alpha = 0.1, correction_method = 'fdr_bh'):
    #     """ Tests the generic function of scimo """
    #     # fimo_file: str, cluster_file: str, cluster_id: str, cluster_gene_id: str,
    #     #                  cluster_p: str, cluster_logfc: str, output_dir: str
    #     scimo = SciMotf('/Users/ariane/Documents/code/sircle/data/S050_CCRCC_Clark_Cell2019/sircle/SIRCLE_A-RNA-protein_TF.csv',
    #                     '/Users/ariane/Documents/code/sircle/data/S050_CCRCC_Clark_Cell2019/sircle/TvN_vis_df.csv', "Regulation_Grouping_2",
    #                     "external_gene_name", "padj_protein", "logFC_protein", self.data_dir,
    #                     fimo_pcol='protein_TF_padj',
    #                     tf_in_dataset=False)
    #     m_df = scimo.run()
    #     # Now we want to assert some simple things since we have validated the output
    #     m_df.to_csv('/Users/ariane/Documents/code/sircle/data/S050_CCRCC_Clark_Cell2019/sircle/SIRCLE_A-RNA-protein_moTF.csv')
    #
    # def test_scimotf2(self):
    #     self.setup_class()
    #     """ Tests the generic function of scimo """
    #     # fimo_file: str, cluster_file: str, cluster_id: str, cluster_gene_id: str,
    #     #                  cluster_p: str, cluster_logfc: str, output_dir: str
    #     fimo_dir = '/Users/ariane/Documents/code/cptac_sircle_admin/data/S050_CCRCC_Clark_Cell2019_SIRCLE/TvN_p005_M10_R05_P02/fimo_100/'
    #
    #     scimo = SciMotf(os.path.join(fimo_dir, 'fimo.tsv'),
    #                     '/Users/ariane/Documents/code/cptac_sircle_admin/data/S050_CCRCC_Clark_Cell2019_SIRCLE/TvN_p005_M10_R05_P02/RCM_TvN-GENES.csv',
    #                     'RegulatoryLabels',
    #                     "external_gene_name",
    #                     "padj_rna",
    #                     "logFC_rna",
    #                     ".",
    #                     cluster_pcutoff=0.1,
    #                     tf_in_dataset=False,
    #                     fimo_pcol="q-value")
    #     m_df = scimo.run()

    def test_doro(self):
        self.setup_class()
        """ Tests the generic function of scimo """
        # fimo_file: str, cluster_file: str, cluster_id: str, cluster_gene_id: str,
        #                  cluster_p: str, cluster_logfc: str, output_dir: str
        scimo = SciMotf(self.fimo, self.cluster, "cluster", "gene_name", "padj", "log2FC", self.data_dir,
                        tf_in_dataset=False)
        m_df = scimo.run()
        # Add dorethea stuff
        dm_df = scimo.add_tf_predictions(m_df, f'{self.data_dir}dorethea_ex.csv', 'TF', 'meanTFChange', "dorethea")
        # Now we want to check that the values are correct
        motif_vals = dm_df['dorethea_meanTFChange'].values
        for i, m in enumerate(dm_df['motif'].values):
            if 'MAZ' in m:
                assert motif_vals[i] == 12.2
            elif 'SP5' in m:
                assert motif_vals[i] == 0.9
            elif 'WT1' in m:
                assert motif_vals[i] == -1.5
