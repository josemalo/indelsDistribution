"""
indelsConsolidation.py
"""

from __future__ import print_function

___author___ = 'Jose Malagon-Lopez'

import sys
import os
import argparse
import collections
import pyfaidx
import HTSeq

# Utilities #
"""
Compute the reverse complement of a genome sequence.
"""
def reverseComplement(seq):
    compl = dict({'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n', '.': '.', '-': '-', '_': '_'})
    out_list = [compl[bp] for bp in seq]
    return ''.join(out_list[::-1])

"""
Read total_counting
"""
def totalCounting(total_count):
    total_counting = {}
    tc = open(total_count, 'Ur')
    lines = tc.readlines()[1:]
    tc.close()

    for line in lines:
        target_name, total_reads, total_ins, total_del, total_mix = line.strip().split('\t')
        total_counting[target_name] = {'total_reads': total_reads, 'total_ins': total_ins, 'total_del': total_del, 'total_mix': total_mix}

    key_list = total_counting.keys()
    key_list.sort()

    return total_counting, key_list

"""
Read mix_out
"""
def MIXout(mix_indels, key_list):
    mix_out = {}
    for key in key_list:
        mix_out[key] = list()

    mi = open(mix_indels, 'Ur')
    items = mi.readlines()[1:]
    mi.close()

    for item in items:
        target_name,target_chromosome,target_start,target_end,target_strand,cleavage_position,target_sequence,total_del,pos_del,len_del,total_ins,pos_ins,len_ins,insertions,indel_seq = item.strip().split('\t')
        mix_outline = [target_name,target_chromosome,target_start,target_end,target_strand,cleavage_position,target_sequence,total_del,pos_del,len_del,total_ins,pos_ins,len_ins,insertions,indel_seq]
        mix_out[target_name].append(mix_outline)
    return mix_out

"""
Read ins_out or del_out
"""
def INDELout(indel_set, key_list):
    set_out = {}
    for key in key_list:
        set_out[key] = list()

    mi = open(indel_set, 'Ur')
    items = mi.readlines()[1:]
    mi.close()

    for item in items:
        target_name,target_chromosome,target_start,target_end,target_strand,cleavage_position,target_sequence,indel_class,indel_start,indel_end,indel_length,insertion,read_name = item.strip().split('\t')
        set_outline = [target_name,target_chromosome,target_start,target_end,target_strand,cleavage_position,target_sequence,indel_class,indel_start,indel_end,indel_length,insertion,read_name]
        set_out[target_name].append(set_outline)
    return set_out

"""
RELATIVE Indel start/end position by strand and w.r.t. the expected cleavage site.
The cleavage site is set as the position in the segment containing the PAM sequence that is closer to the expected cut.
Positive and starting at 1 for every position on the segment containing the PAM sequence.
Negative and starting at -1 for every position on the segment not containing the PAM sequence.
Unlike coordinates which are in the 0-base system, start and end positions are inclusive.
Only one indel at the time.
"""
def indelPosition(indel_pos, strand, indel_start, indel_end, cleavage_pos):
    if strand == '+':
        if indel_pos == 'start_pos':
            indel_rel = int(cleavage_pos) - int(indel_start)
        else:
            indel_rel = int(cleavage_pos) - (int(indel_end) - 1)

        if indel_rel <= 0:
            return abs(indel_rel) + 1
        else:
            return - indel_rel
    else:
        if indel_pos == 'start_pos':
            indel_rel = int(cleavage_pos) - (int(indel_end) - 1)
        else:
            indel_rel = int(cleavage_pos) - int(indel_start)

        if indel_rel >= 0:
            return indel_rel + 1
        else:
            return indel_rel

"""
Remove deletions from Reference Sequence, where an '-' will be store in the deleted position.
"""
def delSequenceMix(genome, chromo, start, end, strand, del_start, del_end):
    pad = 30
    window_start = (int(start) - pad)
    window_end = (int(end) + pad)

    seq = genome[chromo][window_start:window_end].seq.encode().upper()

    if del_start:
        zip_del = zip(del_start, del_end)
        zip_del.sort(key=lambda t: t[0])

        # single block
        if len(del_start) == 1:
            bs = max(int(zip_del[0][0]) - window_start, 0)
            be = min(int(zip_del[0][1]) - window_start, window_end - window_start)
            out_seq = ''.join([seq[:bs], (be - bs) * '-', seq[be:]])

        # multiple blocks: current_index is to be set at the end of each deletion block as the end of the block
        else:
            # first block
            bs = max(int(zip_del[0][0]) - window_start, 0)
            current_index = int(zip_del[0][1]) - window_start
            out_seq = ''.join([seq[:bs], (current_index - bs) * '-'])

            # remaining blocks
            for ds, de in zip_del[1:]:
                gap = int(ds) - window_start - current_index
                block = min((int(de) - int(ds)), (window_end - int(ds)))
                out_seq = out_seq + seq[current_index:(current_index + gap)] + (block * '-')
                current_index += gap + block

            # if necessary, add last block
            if window_end > int(zip_del[-1][1]):
                out_seq += seq[current_index:]
    else:
        out_seq = seq

    if strand == '+':
        return out_seq
    else:
        return reverseComplement(out_seq)


"""
Filter-out indels that do not overlap w. the target sequence.
Filter-out indels of length 1 that are more than two bp away from the expected cleavage site.
For every target, count total number of indels per class.
Consolidate indels as unique per starts/end positions and insertions, if any.
"""
def SingleIndelConsolidation(indel_set, summ_count, genome, variant, replicate, indel_type, output_folder):
    non_allocated, indels_out, total_indels = {}, {}, {}
    #  non_allocated[tag], indels_out[target_name][tag], total_indels[target_name]
    for target_name in summ_count.keys():
        indels_out[target_name], total_indels[target_name] = {}, 0

    for d in indel_set:
        indels_list = indel_set[d]
        for indels in indels_list:
            target_name, target_chromosome, target_start, target_end, target_strand, cleavage_position, target_sequence, indel_class, indel_start, indel_end, indel_length, insertion, read_name = indels
            tag = '_'.join([target_name, str(indel_start), str(indel_end), insertion])

            if int(indel_length) > 1:
                # indel overlap with target sequence
                if (int(target_start) <= int(indel_start) < int(target_end)) or (int(target_start) < int(indel_end) <= int(target_end)) or (int(indel_start) <= int(target_start) and int(target_end) <= int(indel_end)):

                    if tag not in indels_out[target_name]:
                        indel_start_pos = indelPosition('start_pos', target_strand, indel_start, indel_end, cleavage_position)
                        indel_end_pos = indelPosition('end_pos', target_strand, indel_start, indel_end, cleavage_position)

                        if indel_class == 'DEL':
                            del_seq = delSequenceMix(genome, target_chromosome, target_start, target_end, target_strand, [int(indel_start)], [int(indel_end)])
                        else:
                            del_seq = delSequenceMix(genome, target_chromosome, target_start, target_end, target_strand, [], [])

                        indels_out[target_name][tag] = {'target_chromosome': target_chromosome, 'target_start': target_start, 'target_end': target_end,
                                                        'target_strand': target_strand, 'cleavage_position': cleavage_position, 'target_name': target_name,
                                                        'target_sequence': target_sequence, 'indel_class': indel_class, 'indel_start': indel_start, 'indel_end': indel_end,
                                                        'indel_length': indel_length, 'insertion': insertion, 'tag_total': 1,
                                                        'indel_start_pos': indel_start_pos, 'indel_end_pos': indel_end_pos, 'indel_seq': del_seq}
                    else:
                        indels_out[target_name][tag]['tag_total'] += 1

                    if target_name not in total_indels:
                        total_indels[target_name] = 1
                    else:
                        total_indels[target_name] += 1

                # indel does not overlap with target sequence
                else:
                    if tag not in non_allocated:
                        non_allocated[tag] = {'target_chromosome': target_chromosome, 'target_start': target_start, 'target_end': target_end, 'target_strand': target_strand,
                                              'cleavage_position': cleavage_position, 'target_name': target_name, 'target_sequence': target_sequence,
                                              'indel_class': indel_class, 'indel_start': indel_start, 'indel_end': indel_end, 'indel_length': indel_length,
                                              'insertion': insertion, 'tag_total': 1}
                    else:
                        non_allocated[tag]['tag_total'] += 1
            else:
                indel_start_pos = indelPosition('start_pos', target_strand, indel_start, indel_end, cleavage_position)

                # indel of length one must be at most 2 bp from the expected cleavage site
                if abs(indel_start_pos) <= 2:
                    if tag not in indels_out[target_name]:
                        if indel_class == 'DEL':
                            del_seq = delSequenceMix(genome, target_chromosome, target_start, target_end, target_strand, [int(indel_start)], [int(indel_end)])
                        else:
                            del_seq = delSequenceMix(genome, target_chromosome, target_start, target_end, target_strand, [], [])

                        indels_out[target_name][tag] = {'target_chromosome': target_chromosome, 'target_start': target_start, 'target_end': target_end, 'target_strand': target_strand,
                                                        'cleavage_position': cleavage_position, 'target_name': target_name, 'target_sequence': target_sequence,
                                                        'indel_class': indel_class, 'indel_start': indel_start, 'indel_end': indel_end, 'indel_length': indel_length, 'insertion': insertion,
                                                        'tag_total': 1, 'indel_start_pos': indel_start_pos, 'indel_end_pos': (indel_start_pos + 1), 'indel_seq': del_seq}
                    else:
                        indels_out[target_name][tag]['tag_total'] += 1

                    if target_name not in total_indels:
                        total_indels[target_name] = 1
                    else:
                        total_indels[target_name] += 1

                #  indel of length one is more than 2 bp from the expected cleavage site
                else:
                    if tag not in non_allocated:
                        non_allocated[tag] = {'target_chromosome': target_chromosome, 'target_start': target_start, 'target_end': target_end, 'target_strand': target_strand,
                                              'cleavage_position': cleavage_position, 'target_name': target_name, 'target_sequence': target_sequence,
                                              'indel_class': indel_class, 'indel_start': indel_start, 'indel_end': indel_end, 'indel_length': indel_length,
                                              'insertion': insertion, 'tag_total': 1}
                    else:
                        non_allocated[tag]['tag_total'] += 1

    print('Writing Filtered-Out indels.', file=sys.stderr)
    non_allocated_filename = "".join([output_folder, 'NON_allocated/', indel_type, '_', variant, '_', replicate, '_NON_allocated.txt'])

    non_output = os.path.dirname(non_allocated_filename)
    if not os.path.exists(non_output):
        os.makedirs(non_output)

    non_allocated_file = open(non_allocated_filename, 'w')
    print('Variant', 'Replicate', 'Target_Name', 'Target_Chromosome', 'Target_Start', 'Target_End', 'Target_Strand', 'Cleavage_Position',
          'TargetSequence', 'Indel_Class', 'Indel_Start', 'Indel_End', 'Indel_Length', 'Insertion', 'Total',
          file=non_allocated_file, sep='\t')
    non_allocated_file.close()

    with open(non_allocated_filename, 'a') as non_allocated_file:
        for tag in non_allocated:
            outline = variant, replicate, non_allocated[tag]['target_name'], non_allocated[tag]['target_chromosome'], \
                      non_allocated[tag]['target_start'], non_allocated[tag]['target_end'], \
                      non_allocated[tag]['target_strand'], non_allocated[tag]['cleavage_position'], non_allocated[tag]['target_sequence'], \
                      non_allocated[tag]['indel_class'], non_allocated[tag]['indel_start'], non_allocated[tag]['indel_end'], \
                      non_allocated[tag]['indel_length'], non_allocated[tag]['insertion'], non_allocated[tag]['tag_total']
            print(*outline, sep='\t', file=non_allocated_file)

    print('Writing Indel Table.', file=sys.stderr)
    for target_name in indels_out:
        print(target_name, file=sys.stderr)
        indel_filename = "".join([output_folder, target_name, '/', target_name, '_', variant, '_', replicate, '_', indel_type, '.txt'])

        in_output = os.path.dirname(indel_filename)
        if not os.path.exists(in_output):
            os.makedirs(in_output)

        indel_file = open(indel_filename, 'w')
        print('Variant', 'Replicate', 'Target_Name', 'Target_Chromosome', 'Target_Start', 'Target_End', 'Target_Strand', 'Cleavage_Position',
              'TargetSequence', 'Indel_Class', 'Total_Indel_Object', 'PCT', 'Indel_Start', 'Indel_End', 'Indel_Length', 'Insertion',
              'Indel_Start_Position', 'Indel_End_Position', 'Indel_Sequence',
              file=indel_file, sep='\t')
        indel_file.close()

        with open(indel_filename, 'a') as indel_file:
            for tag in indels_out[target_name]:
                outline = variant, replicate, indels_out[target_name][tag]['target_name'], indels_out[target_name][tag]['target_chromosome'], \
                          indels_out[target_name][tag]['target_start'], indels_out[target_name][tag]['target_end'], \
                          indels_out[target_name][tag]['target_strand'], indels_out[target_name][tag]['cleavage_position'], \
                          indels_out[target_name][tag]['target_sequence'], indels_out[target_name][tag]['indel_class'], \
                          indels_out[target_name][tag]['tag_total'], \
                          round(100 * int(indels_out[target_name][tag]['tag_total']) / float(summ_count[target_name]['total_reads']), 3), \
                          indels_out[target_name][tag]['indel_start'], indels_out[target_name][tag]['indel_end'], indels_out[target_name][tag]['indel_length'], \
                          indels_out[target_name][tag]['insertion'], \
                          indels_out[target_name][tag]['indel_start_pos'], indels_out[target_name][tag]['indel_end_pos'], indels_out[target_name][tag]['indel_seq']
                print(*outline, sep='\t', file=indel_file)

    return indels_out, total_indels


def MixedIndelConsolidation(mix_set, summ_count, genome, variant, replicate, output_folder):
    non_allocated, indels_out, total_indels, indels_positions, general_info = {}, {}, {}, {}, {}
    #  non_allocated[tag], indels_out[target_name][tag], total_indels[target_name], indels_positions[target_name][tag], general_info[target_name]

    for target_name in summ_count.keys():
        indels_out[target_name], indels_positions[target_name] = {}, {}

    for d in mix_set:
        indels_list = mix_set[d]
        for item in indels_list:
            target_name,target_chromosome,target_start,target_end,target_strand,cleavage_position,target_sequence,del_total,del_start,del_length,ins_total,ins_start,ins_length,insertion,indel_nbd_seq = item

            if del_start != '.':
                del_start = del_start.split('_')
                del_start = [int(x) for x in del_start]

                del_len = del_length.split('_')
                del_len = [int(x) for x in del_len]
            else:
                del_start, del_len = [], []

            if ins_start != '.':
                ins_start = ins_start.split('_')
                ins_start = [int(x) for x in ins_start]

                ins_len = ins_length.split('_')
                ins_len = [int(x) for x in ins_len]

                insertions = insertion.split('_')
                insertions = [x for x in insertions]
            else:
                ins_start, ins_len, insertions = [], [], []

            #  make tag as a unique identifier of the indel
            if del_start and not ins_start:
                tag = '_'.join([target_name, '_'.join([str(x) for x in del_start]), '_'.join([str(x) for x in del_len])])
            elif not del_start and ins_start:
                tag = '_'.join([target_name, '_'.join(str(x) for x in [ins_start]), '_'.join([str(x) for x in ins_len]), insertion])
            elif del_start and ins_start:
                tag = '_'.join([target_name, '_'.join([str(x) for x in del_start]), '_'.join([str(x) for x in del_len]), '_'.join([str(x) for x in ins_start]), '_'.join([str(x) for x in ins_len]), insertion])

            #  deletions of length one must be at most 2 bp from the expected cleavage site
            del_one_index = [i for i, x in enumerate(del_len) if x == 1]
            if del_one_index:
                index_to_remove = list()
                for i in del_one_index:
                    srp = indelPosition('start_pos', target_strand, int(del_start[i]), int(del_start[i]) + 1, cleavage_position)
                    if abs(srp) > 2:
                        index_to_remove.append(i)
                index_to_remove.reverse()
                for j in index_to_remove:
                    del del_start[j]
                    del del_len[j]

            #  insertions of length one must be at most 2 bp from the expected cleavage site
            ins_one_index = [i for i, x in enumerate(ins_len) if x == 1]
            if ins_one_index:
                index_to_remove = list()
                for i in ins_one_index:
                    srp = indelPosition('start_pos', target_strand, int(ins_start[i]), int(ins_start[i]) + 1, cleavage_position)
                    if abs(srp) > 2:
                        index_to_remove.append(i)
                index_to_remove.reverse()
                for j in index_to_remove:
                    del ins_start[j]
                    del ins_len[j]
                    del insertions[j]

            #  at least one indel must overlap with target sequence
            zip_del = zip(del_start, del_len)
            zip_del.sort(key=lambda t: t[0])
            del_end = [x + y for x, y in zip_del]
            del_start.sort()
            deletion_zip = zip(del_start, del_end)

            zip_ins = zip(ins_start, ins_len)
            zip_ins.sort(key=lambda t: t[0])
            ins_end = [x + y for x, y in zip_ins]
            ins_start.sort()
            insertion_zip = zip(ins_start, ins_end)

            indel_start = del_start + ins_start
            indel_len = del_len + ins_len
            indel_end = [x + y for x, y in zip(indel_start, indel_len)]
            indel_zip = zip(indel_start, indel_end)

            ga_indel = HTSeq.GenomicArray("auto", stranded=False)
            for i,j in indel_zip:
                ga_indel[HTSeq.GenomicInterval("".join(['chr', target_chromosome]), i, j, ".")] = 1
            target_iv = HTSeq.GenomicInterval("".join(['chr', target_chromosome]), int(target_start), int(target_end))

            overlap_indicator = None
            for iv, val in ga_indel[target_iv].steps():
                if val > 0:
                    overlap_indicator = True

            if overlap_indicator:
                if tag not in indels_out[target_name]:
                    indel_seq = delSequenceMix(genome, target_chromosome, target_start, target_end, target_strand, del_start, del_end)

                    indels_out[target_name][tag] = {'target_name': target_name, 'target_chromosome': target_chromosome, 'target_start': target_start,
                                                    'target_end': target_end, 'target_strand': target_strand, 'cleavage_position': cleavage_position,
                                                    'target_sequence': target_sequence, 'del_total': len(del_start),
                                                    'del_start': "_".join([str(x) for x in del_start]),
                                                    'del_length': "_".join([str(x) for x in del_len]),
                                                    'ins_total': len(ins_start),
                                                    'ins_start': "_".join([str(x) for x in ins_start]),
                                                    'ins_length': "_".join([str(x) for x in ins_len]),
                                                    'insertion': "_".join([str(x) for x in insertions]),
                                                    'indel_nbd_seq': indel_nbd_seq, 'tag_total': 1}

                    del_start_pos = [indelPosition('start_pos', target_strand, x, y, cleavage_position) for x, y in deletion_zip]
                    del_end_pos = [indelPosition('end_pos', target_strand, x, y, cleavage_position) for x, y in deletion_zip]
                    ins_start_pos = [indelPosition('start_pos', target_strand, x, y, cleavage_position) for x, y in insertion_zip]
                    ins_end_pos = [indelPosition('end_pos', target_strand, x, y, cleavage_position) for x, y in insertion_zip]

                    indels_positions[target_name][tag] = [del_start_pos, del_end_pos, ins_start_pos, ins_end_pos, indel_seq]

                    general_info[target_name] = [target_chromosome, target_start, target_end, target_strand, cleavage_position, target_sequence]

                else:
                    indels_out[target_name][tag]['tag_total'] += 1

                if target_name not in total_indels:
                    total_indels[target_name] = 1
                else:
                    total_indels[target_name] += 1
            else:
                if tag not in non_allocated:
                    non_allocated[tag] = {'target_name': target_name, 'target_chromosome': target_chromosome, 'target_start': target_start, 'target_end': target_end,
                                          'target_strand': target_strand, 'cleavage_position': cleavage_position, 'target_sequence': target_sequence,
                                          'del_total': len(del_start), 'del_start': "_".join([str(x) for x in del_start]), 'del_length': "_".join([str(x) for x in del_len]),
                                          'ins_total': len(ins_start), 'ins_start': "_".join([str(x) for x in ins_start]), 'ins_length': "_".join([str(x) for x in ins_len]),
                                          'insertion': "_".join([str(x) for x in insertions]), 'indel_nbd_seq': indel_nbd_seq, 'tag_total': 1}
                else:
                    non_allocated[tag]['tag_total'] += 1

    print('Writing Filtered-Out INDELS.', file=sys.stderr)
    non_allocated_filename = "".join([output_folder, 'NON_allocated/', 'INDEL_', variant, '_', replicate, '_NON_allocated.txt'])

    non_output = os.path.dirname(non_allocated_filename)
    if not os.path.exists(non_output):
        os.makedirs(non_output)

    non_allocated_file = open(non_allocated_filename, 'w')
    print('Variant', 'Replicate', 'Target_Name', 'Target_Chromosome', 'Target_Start', 'Target_End',
          'Target_Strand', 'Cleavage_Position', 'TargetSequence', 'Indel_Class',
          'Total_Deletions', 'Deletion_Start', 'Deletion_Length',
          'Total_Insertions', 'Insertion_Start', 'Insertion_Length',
          'Insertion', 'Indel_NBD_Sequence', 'Total_INDEL', file=non_allocated_file, sep='\t')
    non_allocated_file.close()

    with open(non_allocated_filename, 'a') as non_allocated_file:
        for tag in non_allocated:
            outline = variant, replicate, non_allocated[tag]['target_name'], non_allocated[tag]['target_chromosome'], \
                      non_allocated[tag]['target_start'], non_allocated[tag]['target_end'], non_allocated[tag]['target_strand'], \
                      non_allocated[tag]['cleavage_position'], non_allocated[tag]['target_sequence'], 'INDEL', \
                      non_allocated[tag]['del_total'], non_allocated[tag]['del_start'], non_allocated[tag]['del_length'], \
                      non_allocated[tag]['ins_total'], non_allocated[tag]['ins_start'], non_allocated[tag]['ins_length'], \
                      non_allocated[tag]['insertion'], non_allocated[tag]['indel_nbd_seq'], non_allocated[tag]['tag_total']
            print(*outline, sep='\t', file=non_allocated_file)

    print('Writing INDEL Table.', file=sys.stderr)
    for target_name in indels_out:
        print(target_name, file=sys.stderr)

        indel_filename = "".join([output_folder, target_name, '/', target_name, '_', variant, '_', replicate, '_INDEL.txt'])

        in_output = os.path.dirname(indel_filename)
        if not os.path.exists(in_output):
            os.makedirs(in_output)

        indel_file = open(indel_filename, 'w')
        print('Variant', 'Replicate', 'Target_Name', 'Target_Chromosome', 'Target_Start', 'Target_End',
              'Target_Strand', 'Cleavage_Position', 'TargetSequence', 'Indel_Class',
              'Total_Deletions', 'Deletion_Start', 'Deletion_Length',
              'Total_Insertions', 'Insertion_Start', 'Insertion_Length',
              'Insertion', 'Indel_NBD_Sequence', 'Total_INDEL', 'PCT_INDEL', file=indel_file, sep='\t')
        indel_file.close()

        with open(indel_filename, 'a') as indel_file:
            for tag in indels_out[target_name]:
                outline = variant, replicate, indels_out[target_name][tag]['target_name'], indels_out[target_name][tag]['target_chromosome'], \
                          indels_out[target_name][tag]['target_start'], indels_out[target_name][tag]['target_end'], indels_out[target_name][tag]['target_strand'], \
                          indels_out[target_name][tag]['cleavage_position'], indels_out[target_name][tag]['target_sequence'], 'INDEL', \
                          indels_out[target_name][tag]['del_total'], indels_out[target_name][tag]['del_start'], indels_out[target_name][tag]['del_length'], \
                          indels_out[target_name][tag]['ins_total'], indels_out[target_name][tag]['ins_start'], indels_out[target_name][tag]['ins_length'], \
                          indels_out[target_name][tag]['insertion'], indels_out[target_name][tag]['indel_nbd_seq'], indels_out[target_name][tag]['tag_total'], \
                          round(100 * int(indels_out[target_name][tag]['tag_total']) / float(summ_count[target_name]['total_reads']), 3)
                print(*outline, sep='\t', file=indel_file)

    return indels_out, total_indels, indels_positions, general_info


"""
Summary of Count table.
Write Visualization tables.
"""
def indelConsolidation(mix_set, del_set, ins_set, summ_count, genome, variant, replicate, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print('Working with DELETIONS.', file=sys.stderr)
    del_out, total_del = SingleIndelConsolidation(del_set, summ_count, genome, variant, replicate, 'DEL', output_folder)

    print('\nWorking with INSERTIONS.', file=sys.stderr)
    ins_out, total_ins = SingleIndelConsolidation(ins_set, summ_count, genome, variant, replicate, 'INS', output_folder)

    print('\nWorking with MIXED INDELS.', file=sys.stderr)
    indels_out, total_indels, indels_positions, general_info = MixedIndelConsolidation(mix_set, summ_count, genome, variant, replicate, output_folder)

    print('\nWriting Count Summary Table.', file=sys.stderr)
    # updating summary count table
    for target_name in summ_count:
        if target_name in total_del:
            summ_count[target_name]['total_filt_del'] = total_del[target_name]
        else:
            summ_count[target_name]['total_filt_del'] = 0
        if target_name in total_ins:
            summ_count[target_name]['total_filt_ins'] = total_ins[target_name]
        else:
            summ_count[target_name]['total_filt_ins'] = 0
        if target_name in total_indels:
            summ_count[target_name]['total_filt_indels'] = total_indels[target_name]
        else:
            summ_count[target_name]['total_filt_indels'] = 0

    total_count_filename = "".join([output_folder, variant, '_', replicate, '_INDELs_Total_Count.txt'])
    total_count_file = open(total_count_filename, 'w')
    print('Variant', 'Replicate', 'Target_Name', 'Total_Reads', 'Total_Unfiltered_Deletions', 'Total_Unfiltered_Insertions',
          'Total_Filtered_Deletions', 'Total_Filtered_Insertions', 'Total_Filtered_Indels', 'Percentage_Filtered_Indels', file=total_count_file, sep='\t')

    key_summ = summ_count.keys()
    key_summ.sort()

    for target_name in key_summ:
        outline = variant, replicate, target_name, summ_count[target_name]['total_reads'], \
                  summ_count[target_name]['total_del'], summ_count[target_name]['total_ins'],  \
                  summ_count[target_name]['total_filt_del'], summ_count[target_name]['total_filt_ins'], \
                  summ_count[target_name]['total_filt_indels'], 100 * round(float(summ_count[target_name]['total_filt_indels']) / int(summ_count[target_name]['total_reads']), 7)
        print(*outline, sep='\t', file=total_count_file)
    total_count_file.close()

    print('\nWriting Visualization Tables.', file=sys.stderr)
    VIZ_data = {}
    for target_name in general_info:
        VIZ_data[target_name] = collections.defaultdict(list)

    # Run only over the targets with indels to be plot
    window_pad = 30

    for target_name in general_info:
        target_chromosome, target_start, target_end, target_strand, cleavage_position, target_sequence = general_info[target_name]

        # Obtain reference sequence in window
        ref_sequence = genome[target_chromosome][(int(target_start) - window_pad):(int(target_end) + window_pad)].seq.encode()
        if target_strand == '+':
            seq_strand = ref_sequence.upper()
        else:
            seq_strand = reverseComplement(ref_sequence).upper()
        general_info[target_name].append(seq_strand)

    # for target_name in indels_out:
        indel_filename = "".join([output_folder, 'VisualizationTables/', target_name, '_', variant, '_', replicate, '_INDEL.bed'])

        viz_output = os.path.dirname(indel_filename)
        if not os.path.exists(viz_output):
            os.makedirs(viz_output)

        indel_file = open(indel_filename, 'w')
        # target_name, read_count, del_count, ins_count, indel_count, indel_PCT, del_start_pos, del_end_pos, ins_start_pos, ins_end_pos, indel_seq,insertion
        print(target_name, 0, 0, 0, 0, 0, '.', '.', '.', '.', general_info[target_name][-1], '.', file=indel_file, sep='\t')
        indel_file.close()
        VIZ_data[target_name][target_name] = {'name': target_name, 'read_count': 0, 'del_count': 0, 'ins_count': 0, 'indel_count': 0, 'indel_PCT': 0,
                                              'del_start_pos': '.', 'del_end_pos': '.', 'ins_start_pos': '.', 'ins_end_pos': '.',
                                              'del_seq': general_info[target_name][-1], 'insertion': '.'}

        with open(indel_filename, 'a') as indel_file:
            for tag in indels_out[target_name]:
                target_name = indels_out[target_name][tag]['target_name']
                read_count = summ_count[target_name]['total_reads']
                del_count = summ_count[target_name]['total_filt_del']
                ins_count = summ_count[target_name]['total_filt_ins']
                indel_count = indels_out[target_name][tag]['tag_total']
                indel_PCT = round(100 * int(indel_count) / float(read_count), 5)
                if len(indels_positions[target_name][tag][0]) > 0:
                    del_start_pos = "_".join(str(x) for x in indels_positions[target_name][tag][0])
                else:
                    del_start_pos = '.'
                if len(indels_positions[target_name][tag][1]) > 0:
                    del_end_pos = "_".join(str(x) for x in indels_positions[target_name][tag][1])
                else:
                    del_end_pos = '.'
                if len(indels_positions[target_name][tag][2]) > 0:
                    ins_start_pos = "_".join(str(x) for x in indels_positions[target_name][tag][2])
                else:
                    ins_start_pos = '.'
                if len(indels_positions[target_name][tag][3]) > 0:
                    ins_end_pos = "_".join(str(x) for x in indels_positions[target_name][tag][3])
                else:
                    ins_end_pos = '.'
                indel_seq = indels_positions[target_name][tag][4]
                if len(indels_out[target_name][tag]['insertion']) > 0:
                    insertion = indels_out[target_name][tag]['insertion']
                else:
                    insertion = '.'

                outline = "_".join([target_name, del_start_pos, del_end_pos, insertion]), read_count, del_count, ins_count, indel_count, indel_PCT, \
                          del_start_pos, del_end_pos, ins_start_pos, ins_end_pos, indel_seq, insertion
                print(*outline, sep='\t', file=indel_file)

                name, read_count, del_count, ins_count, indel_count, indel_PCT, del_start_pos, del_end_pos, ins_start_pos, ins_end_pos, del_seq, insertion = outline
                VIZ_data[target_name][name] = {'name': name, 'read_count': int(read_count), 'del_count': int(del_count), 'ins_count': int(ins_count),
                                               'indel_count': int(indel_count), 'indel_PCT': indel_PCT, 'del_start_pos': del_start_pos,
                                               'del_end_pos': del_end_pos, 'ins_start_pos': ins_start_pos, 'ins_end_pos': ins_end_pos,
                                               'del_seq': del_seq, 'insertion': insertion}

    return del_out, ins_out, indels_out, summ_count, general_info, VIZ_data

def main():
    parser = argparse.ArgumentParser(description='Consolidate "adequate" indels that overlap the targets')
    parser.add_argument('--mix_indels', help='Mixed indels table from indelsGathering.py', required=True)
    parser.add_argument('--del_indels', help='Deletions table from indelsGathering.py', required=True)
    parser.add_argument('--ins_indels', help='Insertions table from indelsGathering.py', required=True)
    parser.add_argument('--total_count', help='Summary count table from indelsGathering.py', required=True)
    parser.add_argument('--genome_reference', help='Indexed Genome Reference', required=True)
    parser.add_argument('--variant', help='Type of protein used', required=True)
    parser.add_argument('--replicate', help='Replicate number', default=1)
    parser.add_argument('--output_folder', help='Output folder', required=True)
    args = parser.parse_args()

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    genome = pyfaidx.Fasta(args.genome_reference)
    summ_count, key_list = totalCounting(args.total_count)

    mix_set = MIXout(args.mix_indels, key_list)
    del_set = INDELout(args.del_indels, key_list)
    ins_set = INDELout(args.ins_indels, key_list)

    print('*** Running indelsConsolidation ***', file=sys.stderr)
    print(' %s %s' %(args.variant, args.replicate), file=sys.stderr)

    indelConsolidation(mix_set, del_set, ins_set, summ_count, genome, args.variant, args.replicate, args.output_folder)


if __name__ == "__main__":
    main()
