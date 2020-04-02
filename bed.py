import sys

f = sys.argv[1]

with open('refseq_GRCh37_all.bed', 'w+') as outfile:
    with open(f) as infile:
        lines = [l.strip().split('\t') for l in infile.readlines()]
    for line in lines:
        transcript = line[0]
        chrom = line[1]
        start = int(line[3]) - 1000
        end = int(line[4]) + 1000
        outfile.write(
            '{chrom}\t{start}\t{end}\t{name}\n'.format(
                chrom=chrom,
                start=start,
                end=end,
                name=transcript
            )
        )
