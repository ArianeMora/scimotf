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

import argparse
import os
import sys

from sciutil import SciUtil

from scimotf import __version__
from scimotf import SciMotf


def print_help():
    lines = ['-h Print help information.']
    print('\n'.join(lines))


def run(args):
    scimo = SciMotf(args.f, args.c, args.id, args.gene, args.cp, args.logfc, args.o,
                    bg_cluster=args.bg, fimo_pcol=args.fp,
                    fimo_pcutoff=args.fp_c, cluster_pcutoff=args.cp_c, min_genes_in_cluster=args.mg,
                    tf_in_dataset=args.tf, alpha=args.a, correction_method=args.fdr)
    m_df = scimo.run()


def gen_parser():
    parser = argparse.ArgumentParser(description='scimotf')
    parser.add_argument('--f', type=str, help='TSV file output from FIMO')
    parser.add_argument('--c', type=str, help='Output file (csv) from clustering (e.g. from sciRCM)')
    parser.add_argument('--id', type=str, help='Column with cluster annotations in the cluster file.')
    parser.add_argument('--gene', type=str, help='Column name of the gene name from the cluster file.')
    parser.add_argument('--o', type=str, help='Output directory.')

    parser.add_argument('--cp', type=str, help='Column name of the padj value in the cluster file.')
    parser.add_argument('--cp_c', default=1.0, type=float, help='Column name of the logFC value in'
                                                                ' the cluster file (otherwise returns all results).')
    parser.add_argument('--logfc', type=str, help='Column name of the logFC value in the cluster file.')
    parser.add_argument('--bg', type=str, default=None, help='BG cluster row value (e.g. if you have a pre-annotated '
                                                             'cluster name)')
    parser.add_argument('--fp', type=str, default='q-value', help='Column name of the padj value in the FIMO file.')
    parser.add_argument('--fp_c', default=0.05, type=float, help='Cutoff for the pvalue column in the FIMO file (i.e'
                                                                 'the cutoff for the motifs.)')
    parser.add_argument('--mg', type=int, default=3, help='Minimum number of genes per cluster (i.e. too small gene '
                                                          'sets should be removed to improve statistical sig.)')
    parser.add_argument('--tf', type=bool, default=True, help='Whether or not the TF needs to be in the dataset ('
                                                              'recommended as True for RNAseq, and False for protein)')
    parser.add_argument('--a', type=float, default=0.1, help='Alpha value for the FET (i.e looking for abundances in '
                                                             'the clusters).')
    parser.add_argument('--fdr', type=str, default="fdr_bh", help='FDR method for merging the multiple pvalues.')
    return parser


def main(args=None):
    parser = gen_parser()
    u = SciUtil()
    if args:
        sys.argv = args
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f'scimotf v{__version__}')
        sys.exit(0)
    else:
        print(f'scimotf v{__version__}')
        args = parser.parse_args(args)
        # Validate the input arguments.
        if not os.path.isfile(args.f):
            u.err_p([f'The fimo file could not be located, file passed: {args.f}'])
            sys.exit(1)
        if not os.path.isfile(args.c):
            u.err_p([f'The input file that defines the clusters could not be located, file passed: {args.c}'])
            sys.exit(1)

        # Otherwise we have need successful so we can run the program
        u.dp(['Running scimotf: ',
              '\nWith FIMO file: ', args.f,
              '\nAnd cluster file: ', args.c,
              '\nCluster name/id column:', args.id,
              '\nCluster gene ID column: ', args.gene,
              '\nCluster p-adj column:', args.cp,
              '\nCluster logFC: ', args.logfc,
              '\nCluster p-adj cutoff:', args.cp_c,
              '\nBackground cluster: ', args.bg,
              '\nFIMO p-adj column: ', args.fp,
              '\nFIMO p-adj cutoff: ', args.fp_c,
              '\nMinimum number of genes per cluster: ', args.mg,
              '\nWhether of not TF is required in the dataset: ', args.tf,
              '\nFET alpha: ', args.a,
              '\nFET multiple p correction method: ', args.fdr,
              ])
        # RUN!
        run(args)
    # Done - no errors.
    sys.exit(0)


if __name__ == "__main__":
    main()
    # ----------- Example below -----------------------
    """
    root_dir = '../'
    main(["--f", f'{root_dir}tests/data/fimo_sml.tsv',
          "--c", f'{root_dir}tests/data/cluster.csv',
          "--o", f'{root_dir}tests/data/',
          "--id", "cluster",
          "--g", "gene_name",
          "--cp", "padj",
          "--logfc", "log2FC"
          ])
    """
