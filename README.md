# rnaseq
Basic perl script to align RNA-seq single-end fastq files with STAR and then extract featurecounts with featurecounts.
This is designed to be a convenience for RNA-Seq processing, and to save time writing scripts.

## examples
### basic usage
#### specify an annotation file and a STAR genome directory

```shell
perl rnaseq.pl -a ~/GENE_DATA/gencode.v26.annotation.gtf -g ~/GENE_DATA/hg38
```
#### only do debugging, just to see the commands that rnaseq.pl will output

```shell
perl rnaseq.pl -a  ~/GENE_DATA/gencode.v26.annotation.gtf -g ~/GENE_DATA/hg38 -d
```
