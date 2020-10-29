"""
>>> exter = Exter()
>>> exter.add_exon('mygene', [10, 20])
>>> exter.add_exon('mygene', [12, 18])
>>> exter.add_exon('mygene', [30, 40])
>>> exter.add_exon('anothergene', [32, 40])
>>> exter.union()
>>> exter.localize(12, 13)
[{'exon': 0, 'start': 2, 'gene': 'mygene', 'stop': 3, 'strand': '+'}]
>>> exter.localize(32, 33)
[{'exon': 0, 'start': 0, 'gene': 'anothergene', 'stop': 1, 'strand': '+'}, {'exon': 1, 'start': 12, 'gene': 'mygene', 'stop': 13, 'strand': '+'}]
>>> exter.exons('mygene')
[(10, 20), (30, 40)]
>>> exter.globalize('mygene', (9, 11))
[(19, 20), (30, 31)]
>>> exter.globalize('mygene', (30, 41)) # should this raise an error?
[]
>>> exter = Exter()
>>> exter.add_exon('atm', [1000,1100])
>>> exter.add_exon('atm', [5000,5100])
>>> exter.add_exon('atm', [7000,7100])
>>> exter.exons('atm')
[(1000, 1100), (5000, 5100), (7000, 7100)]
>>> exter.localize(5050, 5051)
[{'exon': 1, 'start': 150, 'gene': 'atm', 'stop': 151, 'strand': '+'}]
>>> exter.localize(7050, 7051)
[{'exon': 2, 'start': 250, 'gene': 'atm', 'stop': 251, 'strand': '+'}]
>>> exter.localize(7099, 7100)
[{'exon': 2, 'start': 299, 'gene': 'atm', 'stop': 300, 'strand': '+'}]
>>> exter.localize(1000, 1001)
[{'exon': 0, 'start': 0, 'gene': 'atm', 'stop': 1, 'strand': '+'}]
>>> exter.add_exon('atmm', [1000,1100], reverse_strand=True)
>>> exter.add_exon('atmm', [5000,5100], reverse_strand=True)
>>> exter.localize(5050, 5051)
[{'gene': 'atm', 'start': 150, 'stop': 151, 'strand': '+', 'exon': 1}, {'gene': 'atmm', 'start': 50, 'stop': 51, 'strand': '-', 'exon': 0}]
"""
from __future__ import print_function
from collections import defaultdict
from interlap import InterLap
import gzip
import io
import sys

def read_gff(gff_path, CDS="CDS", union=True):
    """
    read a gff and return a dictionary of chrom -> exter
    """
    D = {'gene': {}, CDS: {}, 'mRNA': {}, 'transcript': {}}
    parents = {}
    fh = io.TextIOWrapper(gzip.open(gff_path, 'r')) if gff_path.endswith(".gz") else open(gff_path, 'r')
    for toks in (x.strip().split("\t") for x in fh if x[0] != '#'):
        if len(toks) > 9: continue
        info =  dict([list(val.split('=')) for val in toks[8].split(';')])
        ftype = toks[2]

        if "ID" in info:
            id = info["ID"].split(":", 1)[-1]
        else:
            id = None

        if 'Parent' in info:
            parent=info['Parent'].split(':', 1)[-1]
        else:
            parent = None
        if id is not None and parent is not None:
            parents[id] = parent

        if not ftype in D:
            continue
        strand = toks[6]
        if ftype == 'gene':
            name = info["Name"]

            if id in D['gene']:
                assert toks[0] == D['gene'][id]["chrom"], toks
                #print(info["Name"], toks, D['gene'][info["Name"]])
                D[ftype][id]['start'] = min(D[ftype][id]['start'], int(toks[3]))
                D[ftype][id]['stop'] = max(D[ftype][id]['stop'], int(toks[4]))

                continue
            D[ftype][id] = dict(chrom=toks[0], start=int(toks[3]),
                    stop=int(toks[4]), strand=strand, name=name)
        else:
            if ftype == 'mRNA': ftype = 'transcript'
            if not id in D[ftype]:
                D[ftype][id] = []
            D[ftype][id].append(dict(chrom=toks[0], start=int(toks[3]),
                stop=int(toks[4]), strand=strand, parent=parent))


    tr = D['transcript']
    result = {}

    for id, exons in D[CDS].items():
        for exon in exons:
            try:
                transr = tr[exon['parent']]
            except:
                # hack it into same structure as above.
                transr = [{"parent": parents[exon['parent']]}]
            try:
                gene = D['gene'][transr[0]['parent']]
            except:
                #print(parents[transr[0]['parent']])
                print("[exter] gene: %s not found" % transr[0]['parent'], file=sys.stderr)
            if not exon['chrom'] in result:
                result[exon['chrom']] = Exter()
            result[exon['chrom']].add_exon(gene['name'], (exon['start'], exon['stop']))

    if union:
        for chrom in result:
            result[chrom].union()

    return result

class Exter(object):
    def __init__(self):
        # TODO use gff
        self._exons = defaultdict(list)
        # this helps find names from positions
        self._posns = InterLap()
        self._dirty = False
        self._reverse = set()

    def add_exon(self, name, start_end, reverse_strand=False):
        self._dirty = True
        self._posns.add((start_end[0], start_end[1], name))
        if reverse_strand:
            if name in self._exons and not name in self._reverse:
                raise Exception(name + " already in as + strand")
            self._reverse.add(name)
        elif name in self._reverse:
            raise Exception(name + " already in as - strand")
        self._exons[name].append((start_end[0], start_end[1]))

    def exons(self, name):
        return self._exons[name]

    def _sort(self):
        for name in self._exons:
            self._exons[name].sort()
        self._dirty = False

    def union(self):
        """
        collapse overlapping exons from the same gene to their union
        """
        self._posns = InterLap()
        for name in self._exons:
            if self._dirty:
                self._exons[name].sort()
            B = self._exons[name]
            before = list(self._exons[name])
            reduced = [self._exons[name][0]]
            for i, b in enumerate(self._exons[name][1:]):
                a = reduced[-1]
                if a[1] >= b[0]:
                    reduced[-1] = (reduced[-1][0], max(reduced[-1][1], b[1]))
                else:
                    reduced.append(b)
            self._exons[name] = reduced
            for r in reduced:
                self._posns.add((r[0], r[1], name))
        self._dirty = False

    def localize(self, start, end):
        if self._dirty: self._sort()

        names = set()
        for p in self._posns.find((start, end)):
            names.add(p[2])

        result = []
        for name in names:
            exons = self._exons[name]
            rev = name in self._reverse
            if rev:
                exons = exons[::-1]

            off = exons[0][0]
            for i, (estart, estop) in enumerate(exons):
                if i > 0:
                    # add last intron
                    #print(off, estart, exons[i-1], estart - exons[i-1][1])
                    off += (estart - exons[i-1][1])

                if estop < start:
                    continue
                if estart > end:
                    continue
                result.append({'gene': name, 'start': start - off, 'stop': end - off,
                    'strand': ('-' if rev else '+'), 'exon': i})

        return result

    def globalize(self, name, start_end):
        """
        """
        if self._dirty: self._sort()
        start, end = start_end
        off = self._exons[name][0][0]

        result = []
        for estart, estop in self._exons[name]:
            if end + off < estart: continue
            if start + off > estop: break
            result.append((max(estart, off + start), min(estop, off + end)))
            off += estop - estart
        return result

    @property
    def exons(self):
        for (s, e, gene) in self._posns:
            yield (s, e, gene)




if __name__ == "__main__":
    import doctest
    print(doctest.testmod())
    d = read_gff(sys.argv[1])
    for chrom in d:
        for (s, e, gene) in d[chrom].exons:
            print("%s\t%d\t%d\t%s" % (chrom, s, e, gene))
