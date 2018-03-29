"""
indelsClosestPrimer.py
"""

from __future__ import print_function

___author___ = 'Jose Malagon-Lopez'

import re
import pyfaidx
import sys
import argparse


# Utilities #
"""
Compute the reverse complement of a genome sequence.
"""
def reverseComplement(seq):
    compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', '.': '.', '-': '-', '_': '_'})
    out_list = [compl[bp] for bp in seq]
    return ''.join(out_list[::-1])

"""
Get smallest distance, if any.
"""
def minDistance(primers_list, pos_index, direction):
    if direction == 'forward':
        dp = [pos_index - x for x in primers_list]
        dp = [x for x in dp if x > 0]
    else:
        dp = [x - pos_index - 1 for x in primers_list]
        dp = [x for x in dp if x > 0]
    if len(dp) > 0:
        return dp[dp.index(min(dp))]
    else:
        return None

# Main Function #

"""
Find closest primer to a given sequence.
"""
def closestPrimer(primer_design, genome, outfile):
    guide_dict = {}
    primers_dict = {}

    pad = 500

    primers = open(primer_design, 'rU')
    for line in primers:
        primer_seq, primer_number, target_name, target_chromosome, target_start, target_end, target_strand, target_sequence = line.strip().split('\t')

        if target_name not in guide_dict:
            if target_strand == '+':
                cleavage_position = int(target_end) - 6
            else:
                cleavage_position = int(target_start) + 5

            guide_dict[target_name] = {'target_chromosome': target_chromosome, 'target_start': int(target_start), 'target_end': int(target_end), 'target_strand': target_strand, 'target_sequence': target_sequence,
                                       'cleavage_position': cleavage_position, 'primer_number': primer_number, 'target_name': target_name}
            primers_dict[target_name] = [primer_seq]
        else:
            primers_dict[target_name].append(primer_seq)
    primers.close()

    for target_name in primers_dict:
        target_chromosome = guide_dict[target_name]['target_chromosome']
        target_start = int(guide_dict[target_name]['target_start'])
        target_end = int(guide_dict[target_name]['target_end'])
        target_sequence = guide_dict[target_name]['target_sequence']

        nbd = genome[target_chromosome][(target_start - pad):(target_end + pad)].seq.encode().upper()
        start_index = pad
        end_index = pad + len(target_sequence)

        p1, p2 = primers_dict[target_name]
        f1, r1, f2, r2 = list(), list(), list(), list()

        mf1 = re.finditer(p1[:20], nbd)
        mr1 = re.finditer(reverseComplement(p1[:20]), nbd)

        mf2 = re.finditer(p2[:20], nbd)
        mr2 = re.finditer(reverseComplement(p2[:20]), nbd)

        # collect positions in nbd sequence of primers occurrence(s)
        for m in mf1:
            f1.append(m.start())
        for m in mr1:
            r1.append(m.start())
        for m in mf2:
            f2.append(m.start())
        for m in mr2:
            r2.append(m.start())

        # determine closets occurrence of primer, either forward or reverse, in nbd
        closest_f1 = minDistance(f1, start_index, 'forward')
        closest_r1 = minDistance(r1, end_index, 'backward')

        #stop in case both forward and reverse primers were not found
        if closest_f1 is None and closest_r1 is None:
            print('Primers for target %s seems to be incorrect.' % target_name, file=sys.stderr)
            sys.exit()

        closest_f2 = minDistance(f2, start_index, 'forward')
        closest_r2 = minDistance(r2, end_index, 'backward')

        #stop in case both forward and reverse primers were not found
        if closest_f2 is None and closest_r2 is None:
            print('Primers for target %s seems to be incorrect.' % target_name, file=sys.stderr)
            sys.exit()

        # choose closest primer occurrence
        if ((closest_f1 and closest_r1 and closest_f2 and closest_r2) is not None and closest_f1 < closest_r1 and closest_f2 > closest_r2) \
                or ((closest_f1 and closest_r1 and closest_r2) is not None and closest_f1 < closest_r1) \
                or ((closest_f1 and closest_f2 and closest_r2) is not None and closest_f2 > closest_r2) \
                or ((closest_f1 and closest_r2) is not None):
            primer_forward = p1
            primer_reverse = p2
            guide_dict[target_name]['PCR_start'] = target_start - closest_f1
            guide_dict[target_name]['PCR_end'] = target_end + closest_r2 + len(primer_reverse) - 1

            if closest_f1 < closest_r2:
                guide_dict[target_name]['primer_strand'] = '+'
                guide_dict[target_name]['primer'] = primer_forward
            else:
                guide_dict[target_name]['primer_strand'] = '-'
                guide_dict[target_name]['primer'] = primer_reverse

        elif ((closest_f1 and closest_r1 and closest_f2 and closest_r2) is not None and closest_f1 > closest_r1 and closest_f2 < closest_r2) \
                or ((closest_f1 and closest_r1 and closest_f2) is not None and closest_f1 > closest_r1) \
                or ((closest_r1 and closest_f2 and closest_r2) is not None and closest_f2 < closest_r2) \
                or ((closest_r1 and closest_f2) is not None):
            primer_forward = p2
            primer_reverse = p1
            guide_dict[target_name]['PCR_start'] = target_start - closest_f2
            guide_dict[target_name]['PCR_end'] = target_end + closest_r1 + len(primer_reverse) - 1

            if closest_r1 < closest_f2:
                guide_dict[target_name]['primer_strand'] = '-'
                guide_dict[target_name]['primer'] = primer_reverse
            else:
                guide_dict[target_name]['primer_strand'] = '+'
                guide_dict[target_name]['primer'] = primer_forward

    output_filename = open(''.join([outfile, '.txt']), 'w')
    print('Closest_Primer_Sequence', 'Primer_Number', 'Target_Name', 'Chromosome',
          'Target_Start', 'Target_End', 'Target_Strand', 'Cleavage_Position',
          'TargetSequence', 'PCR_Start', 'PCR_End', 'Primer_Strand',
          file=output_filename, sep='\t')
    output_filename.close()

    # sort by primer number
    primers = [guide_dict[target_name] for target_name in guide_dict.keys()]
    primers = sorted(primers, key=lambda x: (x['primer_number']))

    with open(''.join([outfile, '.txt']), 'a') as output_filename:
        for primer in primers:
            out = [primer['primer'], primer['primer_number'], primer['target_name'], primer['target_chromosome'],
                   primer['target_start'], primer['target_end'], primer['target_strand'], primer['cleavage_position'],
                   primer['target_sequence'], primer['PCR_start'], primer['PCR_end'], primer['primer_strand']]
            print(*out, sep='\t', file=output_filename)

    with open(''.join([outfile, '.bed']), 'a') as output_filename:
        for primer in primers:
            out = [primer['primer'], primer['target_name'], primer['target_chromosome'], primer['target_start'], primer['target_end'],
                   primer['target_strand'], primer['cleavage_position'], primer['target_sequence']]
            print(*out, sep='\t', file=output_filename)


def main():
    parser = argparse.ArgumentParser(description='Select the closets primer to the expected cleavage site for every target.')
    parser.add_argument('--primer_file', help='Tab separated. Two lines for each target', required=True)
    parser.add_argument('--genome_reference', help='path_to/genome_reference', required=True)
    parser.add_argument('--outfile_basename', help='path_to_output_folder/out_file_basename', required=True)
    args = parser.parse_args()

    genome = pyfaidx.Fasta(args.genome_reference)

    closestPrimer(args.primer_file, genome, args.outfile_basename)

if __name__ == "__main__":
    main()
