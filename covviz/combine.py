import sys
import gzip
import io
import argparse

def gene_regions(bed):

    fh = io.TextIOWrapper(gzip.open(bed, 'r')) if bed.endswith(".gz") else open(bed, 'r')
    for line in fh:
        toks = line.split("\t", 6)
        if len(toks) > 5:
            raise Exception("unknown bed format in %d, expected chrom,start,end,[name,]depth" % bed)
        yield (toks[0], int(toks[1]), int(toks[2]),
                float(toks[4]) if len(toks) == 5 else float(toks[3]),
                toks[3] if len(toks) == 5 else "")
    fh.close()

def combine(beds):
    beds = [gene_regions(b) for b in beds]

    while True:
        try:
            n = next(beds[0])
        except StopIteration:
            for i in range(1, len(beds)):
                try:
                    n = next(beds[i])
                    raise Exception("unequal numbers of lines from file: %d" % i)
                except StopIteration:
                    pass
            return

        row = [n[0], n[1], n[2]]
        if n[4] != "": row.append(n[4])
        row.append(n[3])
        for i in range(1, len(beds)):
            n = next(beds[i])
            assert n[0] == row[0], ("expected identical positions. error in file %d with row %s" % (i, n)) 
            assert n[1] == row[1], ("expected identical positions. error in file %d with row %s" % (i, n)) 
            assert n[2] == row[2], ("expected identical positions. error in file %d with row %s" % (i, n)) 
            row.append(n[3])
        yield row


def combine_main(argv=sys.argv[1:]):
    p = argparse.ArgumentParser("combined bed.gz files")
    p.add_argument("bed", nargs="+", help="coverage files from mosdepth")

    args = p.parse_args(argv)
    print("n bed files: ", len(args.bed), file=sys.stderr)

    for row in combine(args.bed):
        print("%s\t%s\t%s\t%s\t%s" % (row[0], row[1], row[2], row[3],
            "\t".join("%.1f" % d for d in row[4:])))










if __name__ == "__main__":
    combine_main()

