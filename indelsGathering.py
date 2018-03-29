"""
indelsGathering.py
"""

from __future__ import print_function

___author___ = 'Jose Malagon-Lopez'

import sys
import os
import argparse
import pysam
import pyfaidx
import Levenshtein

# Utilities #
"""
Compute the reverse complement of a genome sequence.
"""
def reverseComplement(seq):
    compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', '.': '.', '-': '-', '_': '_'})
    out_list = [compl[bp] for bp in seq]
    return ''.join(out_list[::-1])

"""
Get Levenshtein distance between primer and read-edge, removing up to 5 bp from read edge
"""
def LevenshteinEdge(primer, edge):
    for i in range(6):
        if Levenshtein.distance(primer[i:], edge) <= 5:
            return True
            break

"""
Set indel indicator and retrieve nbd sequence and summary properties.
"""
def gatherIndels(read_seq, cigar_tuples, read_start, read_end, target_chromosome, target_start, target_end, genome):
    # neighboring sequence around the target sequence
    pad = 40
    nbd_start = max(target_start - pad, read_start)
    nbd_end = min(target_end + pad, read_end)
    nbd_seq_index = 0
    nbd_seq = genome[target_chromosome][nbd_start:nbd_end].seq.encode().upper()
    nbd_seq_len = len(nbd_seq)

    current_coordinate = read_start  # zero-base coordinate of the block's initial bp
    rq_seq = read_seq[nbd_start - read_start:]  # observed sequence, to be used in retrieving the insertions
    read_index = 0  # index to be used in tracking the insertions through rq_seq

    total_ins, total_del = 0, 0
    insertion_positions, deletion_positions, insertion_length, deletion_length, insertions = list(), list(), list(), list(), list()
    indel_seq = ''

    for t in cigar_tuples:
        operation = t[0]
        length = t[1]

        # only consider the segments overlapping the nbd sequence
        if current_coordinate + length >= nbd_start and current_coordinate < nbd_end:
            if current_coordinate < nbd_start:  # starts out of the nbd sequence but finishes inside
                distance = current_coordinate + length - nbd_start
                current_coordinate = nbd_start
            else:
                distance = length

            if operation == 0:  # match
                indel_seq += nbd_seq[nbd_seq_index:nbd_seq_index + distance]

                current_coordinate += distance
                read_index += distance
                nbd_seq_index += distance

            elif operation == 1:  # insertion
                total_ins += 1

                indel_seq += nbd_seq[nbd_seq_index:min(nbd_seq_index + distance, nbd_seq_len)]
                insertion = rq_seq[read_index:read_index + distance]

                insertion_positions.append(current_coordinate)
                insertion_length.append(distance)
                insertions.append(insertion)

                read_index += distance

            elif operation == 2:  # deletion
                total_del += 1

                indel_seq += '-' * int(distance)

                deletion_positions.append(current_coordinate)
                deletion_length.append(distance)

                current_coordinate += distance
                nbd_seq_index += distance

    indel_indicator = None
    if total_ins > 0 or total_del > 0:
        indel_indicator = True
    if total_ins == 0:
        insertion_positions, insertion_length, insertions = '.', '.', '.'
    if total_del == 0:
        deletion_positions, deletion_length = '.', '.'

    return indel_indicator, indel_seq, total_del, deletion_positions, deletion_length, total_ins, insertion_positions, insertion_length, insertions


# Main function #
"""
Go through adequate reads in BAM file and obtain explicit description of indels.
"""
def storeIndels(bam_file, primer_file, basename, genome, output_folder):

    print('Reading primers file.', file=sys.stderr)
    total_counting, target_dict = {}, {}

    primers = open(primer_file, 'rU')
    for line in primers:
        primer_seq, target_name, target_chromosome, target_start, target_end, target_strand, cleavage_position, target_sequence = line.strip().split('\t')

        total_counting[target_name] = {'total_reads': 0, 'total_ins': 0, 'total_del': 0, 'total_mix': 0}

        target_dict[target_name] = {'target_name': target_name, 'target_chromosome': target_chromosome, 'target_start': int(target_start), 'target_end': int(target_end),
                                    'target_strand': target_strand, 'target_sequence': target_sequence, 'cleavage_position': int(cleavage_position), 'primer_seq': primer_seq[:20]}
    primers.close()

    # to be use to store alignments w. adequate average score and exactly matching a primer[:20] at one end.
    # At the same time we flagged whether or not the paired alignments are associated to the same target
    read_dict = {}

    ins_out, del_out, mix_out = {}, {}, {}
    for target_name in target_dict:
        ins_out[target_name], del_out[target_name], mix_out[target_name] = list(), list(), list()

    print('Working with the sorted BAM file, one primer_sequence at the time.', file=sys.stderr)
    bamfile = pysam.AlignmentFile(bam_file, 'rb')
    for target_name in target_dict:
        target_start = int(target_dict[target_name]['target_start'])
        target_end = int(target_dict[target_name]['target_end'])
        target_chromosome = target_dict[target_name]['target_chromosome']
        primer_seq = target_dict[target_name]['primer_seq']

        # go over all the reads overlapping with a nbd window of +/-50 bp around the target sequence
        for read in bamfile.fetch(target_chromosome, target_start - 50, target_end + 50):
            read_start = int(read.reference_start)
            read_end = read.reference_end
            read_seq = read.query_sequence
            read_name = read.query_name

            # get first 20-bp sequence in the 5`->3` direction of the query sequence
            if read.is_reverse:
                read_edge = reverseComplement(read_seq[-20:])
            else:
                read_edge = read_seq[:20]

            # only consider the read if:
            # the average quality score of its bases is at least 30;
            # is not supplementary;
            # primer sequence agrees w. the beginning of the query sequence (5'->3');
            # alignment is within 500 bp of target site.
            if sum(read.query_qualities) / len(read.query_qualities) >= 30 and not read.is_supplementary and LevenshteinEdge(primer_seq, read_edge):
                target_strand = target_dict[target_name]['target_strand']
                target_name = target_dict[target_name]['target_name']

                if abs(target_start - read_start) < 500:
                    if read_name not in read_dict:
                        read_dict[read_name] = {'target': None, 'red_flag': None}

                    # if read has been recorded for different primers then activate red_flag, otherwise store its associated primer name
                    if read_dict[read_name]['target'] is not None and read_dict[read_name]['target'] != target_name:
                        read_dict[read_name]['red_flag'] = True
                    else:
                        read_dict[read_name]['target'] = target_name

                    if not read_dict[read_name]['red_flag']:
                        total_counting[target_name]['total_reads'] += 1

                        cigar_tuples = read.cigartuples
                        cigar_string = read.cigarstring

                        if cigar_tuples and cigar_string:
                            if cigar_string.count('D') > 0 and cigar_string.count('I') > 0:
                                total_counting[target_name]['total_mix'] += 1

                            indel_indicator, indel_seq, total_del, deletion_positions, deletion_length, total_ins, insertion_positions, insertion_length, insertions = \
                                gatherIndels(read_seq, cigar_tuples, read_start, read_end, target_chromosome, target_start, target_end, genome)

                            if indel_indicator:
                                if total_del > 0:
                                    total_counting[target_name]['total_del'] += total_del

                                    for i in range(total_del):
                                        outline = [target_name, target_chromosome, target_start, target_end, target_strand,
                                                   target_dict[target_name]['cleavage_position'], target_dict[target_name]['target_sequence'],
                                                   'DEL', deletion_positions[i], int(deletion_positions[i]) + int(deletion_length[i]), deletion_length[i],
                                                   '.', read_name]
                                        del_out[target_name].append(outline)

                                if target_strand == '-':
                                    insertions = [reverseComplement(x) for x in insertions]
                                    indel_seq = reverseComplement(indel_seq)

                                if total_ins > 0:
                                    total_counting[target_name]['total_ins'] += total_ins
                                    for i in range(total_ins):
                                        outline = [target_name, target_chromosome, target_start, target_end, target_strand,
                                                   target_dict[target_name]['cleavage_position'], target_dict[target_name]['target_sequence'],
                                                   'INS', insertion_positions[i], int(insertion_positions[i]) + int(insertion_length[i]), insertion_length[i],
                                                   insertions[i], read_name]
                                        ins_out[target_name].append(outline)

                                mix_outline = [target_name, target_chromosome, target_start, target_end, target_strand,
                                               target_dict[target_name]['cleavage_position'], target_dict[target_name]['target_sequence'],
                                               total_del, '_'.join([str(x) for x in deletion_positions]), '_'.join([str(x) for x in deletion_length]),
                                               total_ins, '_'.join([str(x) for x in insertion_positions]), '_'.join([str(x) for x in insertion_length]),
                                               '_'.join(insertions), indel_seq]
                                mix_out[target_name].append(mix_outline)
    bamfile.close()

    print('Writing INDELS tables.', file=sys.stderr)
    # output files w. indels
    deletion_filename = "".join([output_folder, basename, '_deletions.txt'])
    insertion_filename = "".join([output_folder, basename, '_insertions.txt'])
    indel_filename = "".join([output_folder, basename, '_INDELS.txt'])

    deletion_file = open(deletion_filename, 'w')
    print('Target_Name', 'Target_Chromosome', 'Target_Start', 'Target_End', "Target_Strand", 'Cleavage_Position', 'TargetSequence',
          'Indel_Class', 'Indel_Start', 'Indel_End', 'Indel_Length', 'Insertion', 'Read_Name',
          file=deletion_file, sep='\t')
    deletion_file.close()

    insertion_file = open(insertion_filename, 'w')
    print('Target_Name', 'Target_Chromosome', 'Target_Start', 'Target_End', 'Target_Strand', 'Cleavage_Position', 'TargetSequence',
          'Indel_Class', 'Indel_Start', 'Indel_End', 'Indel_Length', 'Insertion', 'Read_Name',
          file=insertion_file, sep='\t')
    insertion_file.close()

    indel_file = open(indel_filename, 'w')
    print('Target_Name', 'Target_Chromosome', 'Target_Start', 'Target_End', "Target_Strand", 'Cleavage_Position', 'TargetSequence',
          'DEL_Total', 'DEL_Positions', 'DEL_Length', 'INS_Total', 'INS_Positions', 'INS_Length',
          'Insertions', 'Indel_NBD_Sequence',
          file=indel_file, sep='\t')
    indel_file.close()

    for target_name in target_dict:
        with open(insertion_filename, 'a') as insertion_file:
            for insertion in ins_out[target_name]:
                print(*insertion, sep='\t', file=insertion_file)

        with open(deletion_filename, 'a') as deletion_file:
            for deletion in del_out[target_name]:
                print(*deletion, sep='\t', file=deletion_file)

        with open(indel_filename, 'a') as indel_file:
            for indel in mix_out[target_name]:
                print(*indel, sep='\t', file=indel_file)

    # summary count table
    key_list = total_counting.keys()
    key_list.sort()

    summary_filename = "".join([output_folder, basename, '_INDELs_Count_Summary.txt'])
    summary_file = open(summary_filename, 'w')
    print('Target_Name', 'Total_Reads', 'Total_Insertions', 'Total_Deletions', 'Total_Mixed_Indels', file=summary_file, sep='\t')
    summary_file.close()

    with open(summary_filename, 'a') as summary_file:
        for target_name in key_list:
            print(target_name, total_counting[target_name]['total_reads'], total_counting[target_name]['total_ins'], total_counting[target_name]['total_del'], total_counting[target_name]['total_mix'],
                  sep='\t', file=summary_file)

    return total_counting, mix_out, ins_out, del_out


def main():
    parser = argparse.ArgumentParser(description='Gather all the indels for the tag-targeted sites in a 40 bp window.')
    parser.add_argument('--bam_file', help='Sorted bam file with the mapped reads', required=True)
    parser.add_argument('--primer_file', help='Tab separated. A single line per target containing the closest primer to it.', required=True)
    parser.add_argument('--basename', help='basename to be used', required=True)
    parser.add_argument('--genome_reference', help='Indexed Genome Reference', required=True)
    parser.add_argument('--output_folder', help='output folder', required=True)
    args = parser.parse_args()

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    genome = pyfaidx.Fasta(args.genome_reference)

    print('*** Running indelsGathering ***', file=sys.stderr)
    storeIndels(args.bam_file, args.primer_file, args.basename, genome, args.output_folder)

if __name__ == "__main__":
    main()
