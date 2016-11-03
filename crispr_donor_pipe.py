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
import sys
import time
from datetime import datetime
import subprocess
import json

args = docopt(__doc__)


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['primer3'], config_data['master'], config_data['settings']


def populate_seq_dict(id_dict, master, err):
    seq_info = {}
    cur = open(master)
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
                seq_info[nm]['info'] = entry[2:]
                break
    for nm in id_dict:
        if id_dict[nm] == 0:
            err.write(nm + ' not found, skipping!\n')
    return seq_info


def create_seq(nm, info):
        l_input_file = temp_dir + nm + '_LEFT_SEQUENCE.txt'
        r_input_file = temp_dir + nm + '_RIGHT_SEQUENCE.txt'
        left = open(l_input_file, 'w')
        left.write('SEQUENCE_ID=' + nm + '\nSEQUENCE_TEMPLATE=' + info['info'][2] + '\nSEQUENCE_PRIMER_REVCOMP='
                   + info['info'][0] + '\nSEQUENCE_TARGET=37,21\n=')
        left.close()
        right = open(r_input_file, 'w')
        right.write('SEQUENCE_ID=' + nm + '\nSEQUENCE_TEMPLATE=' + info['info'][3] + '\nSEQUENCE_PRIMER='
                   + info['info'][1] + '\nSEQUENCE_TARGET=37,21\n=')
        right.close()
        return l_input_file, r_input_file


def run_primer3(input, output, settings, primer3):
    cmd = primer3 + ' -p3_settings_file=' + settings + ' -output=' + output + ' ' + input
    subprocess.call(cmd, shell=True)


def setup_primer3(seq_dict, primer3, settings, temp_dir):
    for nm in seq_dict:
        (l_input_file, r_input_file) = create_seq(nm, seq_dict[nm])
        l_output_file = temp_dir + nm + '_LEFT_PRIMER3_RESULTS.txt'
        r_output_file = temp_dir + nm + '_RIGHT_PRIMER3_RESULTS.txt'
        run_primer3(l_input_file, l_output_file, settings, primer3)
        run_primer3(r_input_file, r_output_file, settings, primer3)


timestamp = str(int(time.mktime(datetime.now().timetuple())))
warnings = open(timestamp + '_warnings.txt', 'w')
tbl = open(timestamp + '_results.xls', 'w')
temp_dir = timestamp + '_TEMP/'
subprocess.call('mkdir ' + temp_dir, shell=True)

(primer3, master, settings) = parse_config(args['<config>'])
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
setup_primer3(seq_dict, primer3, settings, temp_dir)
