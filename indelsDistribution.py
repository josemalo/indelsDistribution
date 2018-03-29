"""
indelsDistribution.py
"""

from __future__ import print_function

___author___ = 'Jose Malagon-Lopez'

import argparse
import sys
import os
import pysam
import pyfaidx
import Levenshtein
import HTSeq
import svgwrite
import collections
import indelsGathering
import indelsConsolidation
import indelsVIZ

"""
Main Function
"""
def indelsDistribution(genome_reference, primer_file, bam_file, variant, replicate, output_folder):

    if not os.path.exists(os.path.dirname(output_folder)):
        os.makedirs(os.path.dirname(output_folder))

    genome = pyfaidx.Fasta(genome_reference)

    print('\n   *** Gathering Indels *** \n', file=sys.stderr)
    basename = '_'.join([variant, replicate])
    unfiltered_folder = ''.join([output_folder, 'Unfiltered_Indels/'])
    if not os.path.exists(os.path.dirname(unfiltered_folder)):
        os.makedirs(os.path.dirname(unfiltered_folder))

    summ_count, mix_set, ins_set, del_set = indelsGathering.storeIndels(bam_file, primer_file, basename, genome, unfiltered_folder)

    print('\n \n  *** Consolidating Indels *** \n', file=sys.stderr)
    consolidated_folder = ''.join([output_folder, 'Consolidated_Indels/'])
    if not os.path.exists(os.path.dirname(consolidated_folder)):
        os.makedirs(os.path.dirname(consolidated_folder))

    del_out, ins_out, indels_out, summ_count, general_info, VIZ_data = indelsConsolidation.indelConsolidation(mix_set, del_set, ins_set, summ_count, genome, variant, replicate, consolidated_folder)

    print('\n \n *** Writting Visualization Plots ***', file=sys.stderr)
    visualization_folder = ''.join([output_folder, 'Visualization/'])
    if not os.path.exists(os.path.dirname(visualization_folder)):
        os.makedirs(os.path.dirname(visualization_folder))

    for target_name in VIZ_data:
        title = '_'.join([target_name, variant, replicate])
        indelsVIZ.indelsVIZ(VIZ_data[target_name], "".join([visualization_folder, target_name]), target_name, title)

def main():
    parser = argparse.ArgumentParser(description='Consolidate pipeline')
    parser.add_argument('--genome_reference', help='Indexed Genome Reference.', required=True)
    parser.add_argument('--primer_file', help='Tab separated. A single line per target containing the closest primer to it.', required=True)
    parser.add_argument('--bam_file', help='Sorted bam file with the mapped reads.', required=True)
    parser.add_argument('--variant', help='Protein type used', required=True)
    parser.add_argument('--replicate', help='Replicate number', default=1)
    parser.add_argument('--output_folder', help='Main output folder, ending with "/".', required=True)
    args = parser.parse_args()

    indelsDistribution(args.genome_reference, args.primer_file, args.bam_file, args.variant, args.replicate, args.output_folder)

if __name__ == "__main__":
    main()
