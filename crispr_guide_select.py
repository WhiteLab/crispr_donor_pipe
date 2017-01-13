"""Parses pasted output from http://crispr.mit.edu/ and select best guide

Usage: ./crispr_guide_select.py (<flank> <pasted>) [options]

Arguments:
    <flank>    100 bp flanked stop codon reference file
    <pasted>      pasted output from browser

Options:
    -h --help

"""
import re
import sys
from docopt import docopt
args = docopt(__doc__)


def rev_comp(seq):
    code = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    new_seq = ''
    for i in xrange(0, len(seq), 1):
        new_seq += code[seq[i]]
    return new_seq[::-1]


def get_loc(cand, seq):
    m = re.search(cand, seq)
    r = 0
    if m is None:
        rev_cand = rev_comp(cand)
        m = re.search(rev_cand, seq)
        r = 1
    # using 0 to represent spanning stop codon, 1 for 3' in reference to stop, or 2 in 5' reference to stop
    try:
        (start, end) = m.span()
    except:
        return 6, 6, 6
    d = start - 102
    if start <= 100 and end >= 102:
        return 0, 0, r
    elif start > 102:
        return 1, d, r
    else:
        d = 100 - end
        return 2, d, r


def create_oligos(seq):
    root = seq[:-3]
    fwd = 'caccg' + root
    rev = 'aaac' + rev_comp(root) + 'c'
    return fwd, rev

seq_list = []
seq_alias = {}

flanks = open(args['<flank>'])
head = next(flanks)
ind = 0
for flank in flanks:
    info = flank.rstrip('\n').split('\t')
    ids = info[0].split(',')
    syms = info[1].split(',')
    seq_list.append(info[2])
    for nm in ids:
        seq_alias[nm] = ind
    for sym in syms:
        if sym != 'None':
            seq_alias[sym] = ind
    ind += 1
flanks.close()
cur = ''
gRNA = {}
temp = {}
for line in open(args['<pasted>']):
    info = line.rstrip('\n').split('\t')
    if info[0][0:5] != 'Guide':
        cur = info[0]
        gRNA[cur] = {}
        temp = {}
    else:
        (cat, dist, strand) = get_loc(info[2], seq_list[seq_alias[cur]])
        if cat == 6 and dist == 6 and strand == 6:
            sys.stderr.write(cur + '\tCould not find guide with reference. Try using the refseq ID; may have isoforms'
                                   ' with difference stop codons\n')
        if cat not in gRNA[cur]:
            gRNA[cur][cat] = {}
        if dist not in gRNA[cur][cat]:
            gRNA[cur][cat][dist] = {}
        gRNA[cur][cat][dist]['name'] = info[0]
        gRNA[cur][cat][dist]['st'] = strand
        gRNA[cur][cat][dist]['seq'] = info[2]
        gRNA[cur][cat][dist]['score'] = info[1]
# print out prioritized gRNAs
print 'Gene\tCategory\tDistance\tSequence\tSense Oligo\tAntisense Oligo\tScore\tSense\tName'
cat_dict = {0: 'STOP', 1: 'UTR', 2: 'CODON'}
st_dict = {0: '+', 1: '-'}
for gene in gRNA:
    f = 0
    for i in xrange(3):
        if i in gRNA[gene] and f == 0:
            f = 1
            for j in sorted(gRNA[gene][i]):
                cur = gRNA[gene][i][j]
                (forward, reverse) = create_oligos(cur['seq'])
                print '\t'.join((gene, cat_dict[i], str(j), cur['seq'], forward, reverse, cur['score'],
                                 st_dict[cur['st']], cur['name']))
                break
