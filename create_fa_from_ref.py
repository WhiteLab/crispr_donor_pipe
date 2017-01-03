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
        break
    else:
        info = line.rstrip('\n').split('\t')
        ids = info[0].split(',')
        syms = info[1].split(',')

        for i in xrange(len(ids)):
            if ids[i] in id_dict:
                id_dict[ids[i]] = 1
                found += 1
                f = 1
                print '>' + ids[i] + '\n' + info[2]
                break
            elif syms[i] in id_dict:
                id_dict[syms[i]] = 1
                found += 1
                print '>' + syms[i] + '\n' + info[2]
                break
ref.close()

if found != id_len:
    for item in id_dict:
        if id_dict[item] == 0:
            sys.stderr.write('Warning, ' + item + ' not found!\n')
else:
    sys.stderr.write('All items found\n')