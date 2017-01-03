"""Creates fasta file using list of gene symbols or refseq IDs from reference file

Usage: ./create_fa_from_ref (<id_list> <ref>) [options]

Arguments:
    <id_list>    new line-seperated list of gene symbols or refseq IDs
    <ref>      tab-separated reference file with header - ids, symbols, sequence

Options:
    -h --help

"""
from docopt import docopt
import sys

id_dict = {}
args = docopt(__doc__)
id_len = 0
for line in open(args['<id_list>']):
    cur = line.rstrip('\n')
    id_dict[cur] = 0
    id_len += 1

ref = open(args['<ref>'])

head = next(ref)
found = 0
for line in ref:
    if found == id_len:
        ref.close()
    else:
        info = line.rstrip('\n').split('\t')
        ids = info[0].split(',')
        syms = info[1].split(',')
        f = 0
        for nm in ids:
            if nm in id_dict:
                id_dict[nm] = 1
                found += 1
                f = 1
                print nm + '\t' + info[2]
                break
        if f != 0:
            for gene in syms:
                if gene in id_dict:
                    id_dict[gene] = 1
                    found += 1
                    f = 1
                    print gene + '\t' + info[2]
                    break
if found != id_len:
    for item in id_dict:
        if id_dict[item] == 0:
            sys.stderr.write('Warning, ' + item + ' not found!\n')
