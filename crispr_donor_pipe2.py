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


def create_seq(nm, info, max_stop, max_start, side):
        input_file = temp_dir + nm + '_' + side + '_SEQUENCE.txt'
        out_file = open(input_file, 'w')
        seq_len = len(info['seq'][1])
        if side == 'LEFT':
            end = str(seq_len - int(max_stop))
            out_file.write('SEQUENCE_ID=' + nm + '\nSEQUENCE_TEMPLATE=' + info['seq'][0]
                            + '\nSEQUENCE_PRIMER_PAIR_OK_REGION_LIST=' + '0,' + max_start + ',' + end + ','
                           + max_stop + '\n=')
            out_file.close()
        else:
            end = str(seq_len - int(max_start))
            out_file.write('SEQUENCE_ID=' + nm + '\nSEQUENCE_TEMPLATE=' + info['seq'][0]
                           + '\nSEQUENCE_PRIMER_PAIR_OK_REGION_LIST=' + '0,' + max_stop + ',' + end + ','
                           + max_start + '\n=')
            out_file.close()
        return input_file


def parse_results(output, forward, reverse, side, gene):
    f_primer = ''
    r_primer = ''
    attr_dict = {'PRIMER_LEFT_0_PROBLEMS': '', 'PRIMER_LEFT_0_TM': '', 'PRIMER_RIGHT_0_PROBLEMS': '',
                 'PRIMER_RIGHT_0_TM': ''}
    f = 0
    fixed = ''
    for result in open(output):
        cur = result.rstrip('\n').split('=')
        if cur[0] in attr_dict:
            attr_dict[cur[0]] = cur[1]
        if side == 'Left' and cur[0] == 'PRIMER_RIGHT_0_SEQUENCE':
            r_primer = cur[1]
            fixed = r_primer
        elif side == 'Left' and cur[0] == 'PRIMER_LEFT_0_SEQUENCE':
            f_primer = cur[1]
            f = 1
        if side == 'Right' and cur[0] == 'PRIMER_LEFT_0_SEQUENCE':
            f_primer = cur[1]
            fixed = f_primer
        elif side == 'Right' and cur[0] == 'PRIMER_RIGHT_0_SEQUENCE':
            r_primer = cur[1]
            f = 1
    return '\t'.join((gene + '.' + side + '.F', forward + f_primer, attr_dict['PRIMER_LEFT_0_PROBLEMS'],
                      attr_dict['PRIMER_LEFT_0_TM'], gene + '.' + side + '.R', reverse + r_primer,
                      attr_dict['PRIMER_RIGHT_0_PROBLEMS'], attr_dict['PRIMER_RIGHT_0_TM'])), f, fixed


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
        l_input_file = create_seq(nm, seq_dict[nm], max_stop, max_start, 'LEFT')
        r_input_file = create_seq(nm, seq_dict[nm], max_stop, max_start, 'RIGHT')
        l_output_file = temp_dir + nm + '_LEFT_PRIMER3_RESULTS.txt'
        r_output_file = temp_dir + nm + '_RIGHT_PRIMER3_RESULTS.txt'
        gene = seq_dict[nm]['gene']
        run_primer3(l_input_file, l_output_file, Lsettings, primer3)
        run_primer3(r_input_file, r_output_file, Rsettings, primer3)
        # parse results, if primer not found, adjust length and try again
        (left_str, left_flag, left_fixed) = parse_results(l_output_file, lf_gibson, lr_gibson, 'Left', gene)
        if left_flag == 0:
            # cur_gc = calc_gc(left_fixed)
            warn = 'No primer for ' + nm + ' left seq found at ' + max_start + ' left and ' + max_stop \
                   + ' from stop codon\n'

            warnings.write(warn)
        (right_str, right_flag, right_fixed) = parse_results(r_output_file, rf_gibson, rr_gibson, 'Right', gene)
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


header = 'RefSeq ID \tDonor Left join F\tDonor Left join F oligo sequence\tLF_Problems\tLF_TM\tDonor Left join R\t' \
         'Donor Left join R oligo sequence\tLR_Problems\tLR_TM\tDonor Right join F\tDonor Right join F oligo sequence' \
         '\tRF_Problems\tRF_TM\tDonor Right join R\tDonor Right join R oligo sequence\tRR_Problems\tRR_TM\n'
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
