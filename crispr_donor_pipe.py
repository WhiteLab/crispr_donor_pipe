"""Creates primers using refseq ID list and master table

Usage: ./crispr_donor_pipe.py (<config> <list> <table>) [options]

Arguments:
    <config>    json formatted config file with refs and tool locations
    <list>      list of refSeq IDs, one per line
    <table>     Table with refseq ID, gene symbol, starting left and right primer sequences and left and right primer regions

Options:
    -h --help
    --lf  LEFT_FORWARD  left forward sequence to prepend
    --lr  LEFT_REVERSE  left reverse sequence to prepend
    --rf  RIGHT_FORWARD  right forward sequence to prepend
    --rr  RIGHT_REVERSE  right reverse sequence to prepend

"""
from docopt import docopt
import gzip
import sys
import time
from datetime import datetime
import subprocess
import json

args = docopt(__doc__)


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['primer3'], config_data['master'], config_data['Lsettings'], config_data['Rsettings'], \
           config_data['LR_len'], config_data['RF_len']


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


def create_seq(nm, info, LR_len, RF_len):
        l_input_file = temp_dir + nm + '_LEFT_SEQUENCE.txt'
        r_input_file = temp_dir + nm + '_RIGHT_SEQUENCE.txt'
        left = open(l_input_file, 'w')
        lr_prime = info['seq'][0][(int(LR_len) * -1):]
        lr_prime = rev_comp(lr_prime)
        left.write('SEQUENCE_ID=' + nm + '\nSEQUENCE_TEMPLATE=' + info['seq'][0] + '\nSEQUENCE_PRIMER_REVCOMP='
                   + lr_prime + '\nSEQUENCE_TARGET=37,21\n=')
        left.close()
        rf_prime = info['seq'][1][:(int(RF_len))]
        right = open(r_input_file, 'w')
        right.write('SEQUENCE_ID=' + nm + '\nSEQUENCE_TEMPLATE=' + info['seq'][1] + '\nSEQUENCE_PRIMER='
                   + rf_prime + '\nSEQUENCE_TARGET=37,21\n=')
        right.close()
        return l_input_file, r_input_file


#def parse_results(output):



def run_primer3(input, output, settings, primer3):
    cmd = primer3 + ' -p3_settings_file=' + settings + ' -output=' + output + ' ' + input
    subprocess.call(cmd, shell=True)


def setup_primer3(seq_dict, primer3, Lsettings, Rsettings, temp_dir, LR_len, RF_len):
    for nm in seq_dict:
        (l_input_file, r_input_file) = create_seq(nm, seq_dict[nm], LR_len, RF_len)
        l_output_file = temp_dir + nm + '_LEFT_PRIMER3_RESULTS.txt'
        r_output_file = temp_dir + nm + '_RIGHT_PRIMER3_RESULTS.txt'
        run_primer3(l_input_file, l_output_file, Lsettings, primer3)
        run_primer3(r_input_file, r_output_file, Rsettings, primer3)


timestamp = str(int(time.mktime(datetime.now().timetuple())))
warnings = open(timestamp + '_warnings.txt', 'w')
tbl = open(timestamp + '_results.xls', 'w')
temp_dir = timestamp + '_TEMP/'
subprocess.call('mkdir ' + temp_dir, shell=True)

(primer3, master, Lsettings, Rsettings, LR_len, RF_len) = parse_config(args['<config>'])
header = 'RefSeq ID \tDonor Left join F\tDonor Left join F oligo sequence\tDonor Left join R\t' \
         'Donor Left join R oligo sequence\tDonor Right join F\tDonor Right join F oligo sequence\t' \
         'Donor Right join F\tDonor Right join F oligo sequence\n'
tbl.write(header)
id_dict = {}
# set up transcript list
for line in open(args['<list>']):
    line = line.rstrip('\n')
    id_dict[line] = 0
# get relevant seqs from table
seq_dict = populate_seq_dict(id_dict, master, warnings)
setup_primer3(seq_dict, primer3, Lsettings, Rsettings, temp_dir, LR_len, RF_len)
