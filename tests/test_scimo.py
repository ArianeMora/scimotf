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
