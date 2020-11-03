covviz provides several helpers to faciliate viewing coverage of exome data.

## gff-exons

### Extracting exons from GFF

To extract exons from a GFF file

```
covviz gff-exons $gff > exons.bed
```

### Exon coverage

Given the exons.bed as created above, use [mosdepth](https://github.com/brentp/mosdepth) to get the
median coverage for each exon (for each sample).


```
mosdepth --use-median -t 2 -x -f $fasta -n --by exons.bed $sample_name $cram'
```

This will create `$sample_name.regions.bed.gz` for each sample.

### Combining Samples

Mosdepth outputs a single file for each sample. The `combine` command makes these into a single file.
**NOTE**: you must send the `$sample.regions.bed.gz`, **not** the per-base file.

```
covviz combine *.regions.bed.gz > combined.bed
```

### Running covviz

Given a `combined` file a above. Then you can run `covviz` like this:

```
covviz
