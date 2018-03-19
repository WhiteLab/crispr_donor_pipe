"""Creates primers using refseq ID list and master table

Usage: ./crispr_donor_pipe2.py (<config> <list> <max_stop> <max_start>) [options]

Arguments:
    <config>    json formatted config file with refs and tool locations
    <list>      list of refSeq IDs, one per line
    <max_stop>    max num nt from stop codon to search from (i.e., 40)
    <max_start>    max num nt from start to search from (i.e. 380)

Options:
    -h --help

"""
from docopt import docopt
import gzip
import time
import subprocess
import json
import os
import re

args = docopt(__doc__)


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['primer3'], config_data['master'], config_data['Lsettings'], config_data['Rsettings'], \
           config_data['lf_gibson'], config_data['lr_gibson'], config_data['rf_gibson'], config_data['rr_gibson'], \
           config_data['GC_trigger']


def populate_seq_dict(id_dict, master, err):
    seq_info = {}
    cur = gzip.open(master)
    next(cur)
    for entry in cur:
        entry = entry.rstrip('\n').split('\t')
        id_list = entry[0].split(',')
        gene_list = entry[1].split(',')
        for i in xrange(0, len(id_list), 1):
            nm = id_list[i]
            if nm in id_dict:
                id_dict[nm] = 1
                seq_info[nm] = {}
                seq_info[nm]['gene'] = gene_list[i]
                seq_info[nm]['seq'] = entry[2:]
                break
    for nm in id_dict:
        if id_dict[nm] == 0:
            err.write(nm + ' not found, skipping!\n')
    return seq_info


def rev_comp(seq):
    code = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    new_seq = ''
    for i in xrange(0, len(seq), 1):
        new_seq += code[seq[i]]
    return new_seq[::-1]


def create_seq(nm, seq, max_stop, max_start, side):
        input_file = temp_dir + nm + '_' + side + '_SEQUENCE.txt'
        out_file = open(input_file, 'w')
        seq_len = len(seq)
        if side == 'LEFT':

            end = str(seq_len - int(max_stop))
            out_file.write('SEQUENCE_ID=' + nm + '\nSEQUENCE_TEMPLATE=' + seq
                           + '\nSEQUENCE_PRIMER_PAIR_OK_REGION_LIST=' + '0,' + max_start + ',' + end + ','
                           + max_stop + '\n=')
            out_file.close()
        else:
            end = str(seq_len - int(max_start))
            out_file.write('SEQUENCE_ID=' + nm + '\nSEQUENCE_TEMPLATE=' + seq
                           + '\nSEQUENCE_PRIMER_PAIR_OK_REGION_LIST=' + '0,' + max_stop + ',' + end + ','
                           + max_start + '\n=')
            out_file.close()
        return input_file


def process_hits(num_res, side, seq, fh, gene):
    temp = {}
    dist = len(seq)
    best_index = dist
    best_dist = dist

    for result in fh:
        data = result.rstrip('\n').split('=')
        m = re.match('PRIMER_([LEFT|RIGHT])_(\d)_([SEQUENCE|TM])', data[0])
        if m:
            (cur_side, cur_hit, cur_info) = (m.group(1), m.group(2), m.group(3))
            if cur_hit not in temp:
                temp[cur_hit] = {}
            if cur_info == 'TM':
                if cur_side == 'LEFT':
                    temp[cur_hit]['l_tm'] = data[1]
                else:
                    # after getting melting point for last hit, no more info needed from file
                    temp[cur_hit]['r_tm'] = data[1]
                    if int(cur_hit) == num_res - 1:
                        fh.close()
            elif side == 'Left' and cur_side == 'RIGHT':
                r_primer = data[1]
                cur_dist = dist - seq.rfind(rev_comp(r_primer)) + len(r_primer)
                if cur_dist < best_dist:
                    best_index = cur_hit
                    best_dist = cur_dist
                temp[cur_hit]['r_primer'] = r_primer
            elif side == 'Left' and cur_side == 'LEFT':
                f_primer = data[1]
                temp[cur_hit]['f_primer'] = f_primer
            if side == 'Right' and cur_side == 'LEFT':
                f_primer = data[1]
                temp[cur_hit]['f_primer'] = f_primer
                cur_dist = dist - seq.find(rev_comp(f_primer))
                if cur_dist < best_dist:
                    best_index = cur_hit
                    best_dist = cur_dist

            elif side == 'Right' and cur_side == 'RIGHT':
                r_primer = data[1]
                temp[cur_hit]['r_primer'] = r_primer
    warnings.write('Best hit for ' + side + ' for ' + gene + ' was ' + best_index + ' (counting from 0) which was '
                   + str(best_dist) + ' away\n')
    return temp[best_index]['f_primer'], temp[best_index]['r_primer'], temp[best_index]['l_tm'], temp[best_index]['r_tm']


def parse_results(output, forward, reverse, side, gene, seq):
    f_primer = ''
    r_primer = ''
    l_tm = ''
    r_tm = ''
    attr_dict = {'PRIMER_LEFT_EXPLAIN': '',  'PRIMER_RIGHT_EXPLAIN': ''}
    f = 0
    fixed = ''
    fh = open(output)

    for result in fh:
        cur = result.rstrip('\n').split('=')
        if cur[0] == 'PRIMER_PAIR_NUM_RETURNED':
            if cur[1] == 0:
                fh.close()
                break
            else:
                num_res = int(cur[1])
                f = 1
                (f_primer, r_primer, l_tm, r_tm) = process_hits(num_res, side, seq, fh, gene)
        if cur[0] in attr_dict:
            attr_dict[cur[0]] = cur[1]

    return '\t'.join((gene + '.' + side + '.F', forward + f_primer, attr_dict['PRIMER_LEFT_EXPLAIN'],
                      l_tm, gene + '.' + side + '.R', reverse + r_primer,
                      attr_dict['PRIMER_RIGHT_EXPLAIN'], r_tm)), f, fixed


def run_primer3(input, output, settings, primer3):
    cmd = '\"' + primer3 + '\" -p3_settings_file=\"' + settings + '\" -output=' + output + ' ' + input
    subprocess.call(cmd, shell=True)


def calc_gc(seq):
    slen = len(seq)
    gc_ct = 0.0
    gc_dict = {'G': 0.0, 'C': 0.0}
    for nt in seq:
        if nt in gc_dict:
            gc_ct += 1.0
            gc_dict[nt] += 1.0
    return gc_ct/slen


def setup_primer3(seq_dict, primer3, Lsettings, Rsettings, temp_dir, max_stop, max_start, lf_gibson, lr_gibson,
                  rf_gibson, rr_gibson, tbl, warnings):
    for nm in seq_dict:
        # L and R seq in array 0 and 1
        l_input_file = create_seq(nm, seq_dict[nm]['seq'][0], max_stop, max_start, 'LEFT')
        r_input_file = create_seq(nm, seq_dict[nm]['seq'][1], max_stop, max_start, 'RIGHT')
        l_output_file = temp_dir + nm + '_LEFT_PRIMER3_RESULTS.txt'
        r_output_file = temp_dir + nm + '_RIGHT_PRIMER3_RESULTS.txt'
        gene = seq_dict[nm]['gene']
        run_primer3(l_input_file, l_output_file, Lsettings, primer3)
        run_primer3(r_input_file, r_output_file, Rsettings, primer3)
        # parse results, if primer not found, adjust length and try again
        (left_str, left_flag, left_fixed) = parse_results(l_output_file, lf_gibson, lr_gibson, 'Left', gene,
                                                          seq_dict[nm]['seq'][0])
        if left_flag == 0:
            # cur_gc = calc_gc(left_fixed)
            warn = 'No primer for ' + nm + ' left seq found at ' + max_start + ' left and ' + max_stop \
                   + ' from stop codon\n'

            warnings.write(warn)
        (right_str, right_flag, right_fixed) = parse_results(r_output_file, rf_gibson, rr_gibson, 'Right', gene,
                                                             seq_dict[nm]['seq'][0])
        if right_flag == 0:
            warn = 'No primer for ' + nm + ' right seq found at ' + max_start + ' from stop codon and ' + max_stop \
                   + ' from stop right\n'

            warnings.write(warn)
        tbl.write(nm + '\t' + left_str + '\t' + right_str + '\n')


# ACTUAL START OF PROGRAM
(primer3, master, Lsettings, Rsettings, lf_gibson, lr_gibson, rf_gibson, rr_gibson, gc_trigger) = \
    parse_config(args['<config>'])
(max_stop, max_start) = (args['<max_stop>'], args['<max_start>'])
timestamp = time.strftime("%Y-%m-%d_%H%M") + '_' + max_stop + '_' + max_start
warnings = open(timestamp + '_warnings.txt', 'w')
tbl = open(timestamp + '_results.xls', 'w')
temp_dir = timestamp + '_TEMP/'
os.mkdir(temp_dir)


header = 'RefSeq ID \tDonor Left join F\tDonor Left join F oligo sequence\tLF_Notes' \
         '\tLF_TM\tDonor Left join R\tDonor Left join R oligo sequence\tLR_Notes\tLR_TM\tDonor Right join F' \
         '\tDonor Right join F oligo sequence\tRF_Notes\tRF_TM\tDonor Right join R\tDonor Right join R oligo sequence' \
         '\tRR_Notes\tRR_TM\n'
tbl.write(header)
id_dict = {}
# set up transcript list
for line in open(args['<list>']):
    line = line.rstrip('\n')
    id_dict[line] = 0
# get relevant seqs from table
seq_dict = populate_seq_dict(id_dict, master, warnings)
setup_primer3(seq_dict, primer3, Lsettings, Rsettings, temp_dir, max_stop, max_start, lf_gibson, lr_gibson, rf_gibson,
              rr_gibson, tbl, warnings)
tbl.close()
warnings.close()
