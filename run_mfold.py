"""Helper wrapper script for crispr_donor_pipe2.py for secondary structure predictions

Usage: ./run_mfold.py (<config> <seg_name> <seq> <tm> <out>) [options]

Arguments:
    <config>    json formatted config file with refs and tool locations
    <seg_name>      short, but descriptive name for sequence
    <seq>    sequence to fold
    <tm>    melting temp to test fold at AS INTEGER, rounded down
    <out>   out dir location, i.e. 2018-03-20_1404_40_360_TEMP/FOLDS

Options:
    -h --help

"""

import subprocess
import json
import os
import sys


def parse_config(config_file):
    config_data = json.loads(open(config_file, 'r').read())
    return config_data['mfold']


def create_fa(out_dir, seg_name, seq):
    if not os.path.isdir(out_dir):
        mk_out = 'mkdir -p ' + out_dir
        subprocess.call(mk_out, shell=True)
    fn = out_dir + seg_name + '.fa'
    fa = open(seg_name + '.fa', 'w')
    fa.write('>' + seg_name + '\n' + seq + '\n')
    fa.close()
    return fn


def run_mfold(config, out_dir, seg_name, seq, tm):
    mfold = parse_config(config)
    fn = create_fa(out_dir, seg_name, seq)
    run_mfold_cmd = mfold + ' SEQ=' + fn + ' NA=DNA T=' + tm
    sys.stderr.write('mfold command ' + run_mfold_cmd + '\n')
    check = subprocess.call(run_mfold_cmd, shell=True)
    if check != 0:
        sys.stderr.write('mfold exited with errors for ' + seg_name + '\n')
        return 1
    else:
        return 0


# The fun starts here
if __name__ == '__main__':
    from docopt import docopt
    args = docopt(__doc__)
    run_mfold(args['<config>'], args['<out>'], args['<seg_name>'], args['<seq>'], args['<tm>'])

