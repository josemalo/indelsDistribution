""" indelsVIZ.py """

from __future__ import print_function

___author___ = 'Jose Malagon-Lopez'

import collections
import svgwrite
import sys
import os
import argparse

"""
Read VIZdata
"""
def VIZdata(infile):
    indels_dict = collections.defaultdict(list)

    f = open(infile, 'rU')
    lines = f.readlines()
    f.close()

    print('Reading File', file=sys.stderr)
    for line in lines:
        name, read_count, del_count, ins_count, indel_count, indel_PCT, del_start_pos, del_end_pos, ins_start_pos, ins_end_pos, del_seq, insertion = line.strip().split('\t')

        indels_dict[name] = {'name': name, 'read_count': int(read_count), 'del_count': int(del_count),
                             'ins_count': int(ins_count),
                             'indel_count': int(indel_count), 'indel_PCT': indel_PCT,
                             'del_start_pos': del_start_pos, 'del_end_pos': del_end_pos,
                             'ins_start_pos': ins_start_pos, 'ins_end_pos': ins_end_pos,
                             'del_seq': del_seq, 'insertion': insertion}
    return indels_dict

def indelsVIZ(indels_dict, outfile, target_name, title):
    print('\n', file=sys.stderr)
    print(target_name, file=sys.stderr)
    # make a list out of them
    indels = [indels_dict[index] for index in indels_dict.keys()]
    indels = sorted(indels, key=lambda x: (x['indel_count']), reverse=True)

    print('Initiating Canvas', file=sys.stderr)
    box_size = 15
    dwg = svgwrite.Drawing(outfile + '.svg', profile='full', size=(u'131%', 100 + len(indels) * (box_size + 1)))

    if title is not None:
        # Define top and left margins
        x_offset = 10
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 30
        y_offset = 20

    # Basic objects
    ref_seq = indels_dict[target_name]['del_seq']

    # Restrict indels_dict to indels only
    indels.remove(indels_dict[target_name])

    # Draw ticks
    upper_lim = len(ref_seq)

    tick_locations = [1, upper_lim]
    tick_locations += range(upper_lim)[::10][1:]

    for x in tick_locations:
        if x <= 46:
            t = x - 48
            dwg.add(dwg.text(str(t), insert=(x_offset + (x - 1.5) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))
        else:
            t = x - 47
            dwg.add(dwg.text(str(t), insert=(x_offset + (x - 1.1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

    print('Drawing Indels')
    # Draw reference sequence row: del_seq, indel_count, indel_PCT, del_length, insertion

    for i, bp in enumerate(ref_seq):
        y = y_offset
        x = x_offset + i * box_size

        if i <= 29 or i >= 53:
            dwg.add(dwg.text(bp, insert=(x, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))

        elif 30 <= i <= 46:
            dwg.add(dwg.text(bp, insert=(x, 2 * box_size + y - 3), fill='red', style="font-size:15px; font-family:Courier", font_weight='bold'))

        else:
            dwg.add(dwg.text(bp, insert=(x, 2 * box_size + y - 3), fill='navy', style="font-size:15px; font-family:Courier", font_weight='bold'))

    dwg.add(dwg.text('Count', insert=(x_offset + box_size * upper_lim + 16, y_offset + box_size - 5),
                     style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('PCT', insert=(x_offset + box_size * upper_lim + 100, y_offset + box_size - 5),
                     style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('DEL_Length', insert=(x_offset + box_size * upper_lim + 165, y_offset + box_size - 5),
                     style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('INS', insert=(x_offset + box_size * upper_lim + 285, y_offset + box_size - 5),
                     style="font-size:15px; font-family:Courier"))

    # Draw indels
    y_offset += 25  # leave some extra space after the reference row

    for j, seq in enumerate(indels):
        del_seq = indels[j]['del_seq']
        k = 0
        y = y_offset + j * box_size

        # draw positions where insertions take place
        if indels[j]['ins_end_pos'] != '.':
            ins = [int(x) for x in indels[j]['ins_start_pos'].split('_')]
            ine = [int(x) for x in indels[j]['ins_end_pos'].split('_')]

            for s, e in zip(ins, ine):
                # set positions
                if s < 0:
                    # set insertion-box length
                    if e > 0:
                        ins_len = abs(s) + e
                    else:
                        ins_len = abs(s) - abs(e) + 1

                    # draw grey box
                    dwg.add(dwg.rect((x_offset + (s + 47) * box_size - 3, box_size + y), (ins_len * box_size, box_size), fill='lightsteelblue'))
                    # draw vertical line
                    dwg.add(dwg.text(u"\u007C", insert=(x_offset + (s + 46) * box_size + 10, y_offset + box_size * (j + 2) - 2),
                                     fill='black', style="font-size:17px; font-family:Impact", font_weight='bolder'))

                else:
                    # set insertion-box length; insertion may go beyond the shown sequence
                    ins_len = min(abs(e) - abs(s) + 1, 38 - abs(s) + 1)

                    # draw grey box
                    dwg.add(dwg.rect((x_offset + (s + 46) * box_size - 3, box_size + y), (ins_len * box_size, box_size), fill='lightsteelblue'))
                    # draw vertical line
                    dwg.add(dwg.text(u"\u007C", insert=(x_offset + (s + 45) * box_size + 10, y_offset + box_size * (j + 2) - 2),
                                     fill='black', style="font-size:17px; font-family:Impact", font_weight='bolder'))

        # draw deletion sequence
        for i, bp in enumerate(del_seq):
            x = x_offset + k * box_size
            if bp == '-':
                dwg.add(dwg.text(u"\u2013", insert=(x, 2 * box_size + y - 4), fill='black', style="font-size:103 px; font-family:Courier", font_weight='bolder'))
                k += 1

            elif i <= 29 or i >= 53:
                dwg.add(dwg.text(bp, insert=(x, 2 * box_size + y - 3), fill='grey', style="font-size:15px; font-family:Courier"))
                k += 1

            elif 30 <= i <= 46:
                dwg.add(dwg.text(bp, insert=(x, 2 * box_size + y - 3), fill='red', style="font-size:15px; font-family:Courier", font_weight='bold'))
                k += 1
            else:
                dwg.add(dwg.text(bp, insert=(x, 2 * box_size + y - 3), fill='navy', style="font-size:15px; font-family:Courier", font_weight='bold'))
                k += 1

        # set deletion(s) length
        del_length = 0
        del_len = []
        if indels[j]['del_end_pos'] != '.':
            ds = [int(x) for x in indels[j]['del_start_pos'].split('_')]
            de = [int(x) for x in indels[j]['del_end_pos'].split('_')]
            for s,e in zip(ds, de):
                if s > 0 and e > 0:
                    del_len.append(str(abs(e) - abs(s) + 1))
                elif s < 0 and e < 0:
                    del_len.append(str(abs(s) - abs(e) + 1))
                else:
                    del_len.append(str(abs(s) + e))
            del_length = '_'.join(del_len)

        indel_count = seq['indel_count']
        reads_text1 = dwg.text(indel_count, insert=(box_size * (upper_lim + 1) + 31, y_offset + box_size * (j + 2) - 2),
                               fill='black', style="font-size:15px; font-family:Courier")
        dwg.add(reads_text1)

        indel_PCT = seq['indel_PCT']
        reads_text2 = dwg.text(indel_PCT, insert=(box_size * (upper_lim + 1) + 91, y_offset + box_size * (j + 2) - 2),
                               fill='black', style="font-size:15px; font-family:Courier")
        dwg.add(reads_text2)

        del_length = del_length
        reads_text3 = dwg.text(del_length, insert=(box_size * (upper_lim + 1) + 193, y_offset + box_size * (j + 2) - 2),
                               fill='black', style="font-size:15px; font-family:Courier")
        dwg.add(reads_text3)

        insertion = seq['insertion']
        reads_text4 = dwg.text(insertion, insert=(box_size * (upper_lim + 1) + 289, y_offset + box_size * (j + 2) - 2),
                               fill='black', style="font-size:15px; font-family:Courier")
        dwg.add(reads_text4)
    dwg.save()

def main():
    parser = argparse.ArgumentParser(description='Visualization plots with multiple field values')
    parser.add_argument('--infile', help='Input file', required=True)
    parser.add_argument('--target_name', help="guide's name", required=True)
    parser.add_argument('--output_folder', help='Output folder', required=True)
    parser.add_argument('--title', help='title plot', required=True)
    args = parser.parse_args()

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    indels_dict = VIZdata(args.infile)

    indelsVIZ(indels_dict, "".join([args.output_folder, args.target_name]), args.target_name, args.title)


if __name__ == "__main__":
    main()
