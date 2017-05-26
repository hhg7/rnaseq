# rnaseq
Basic perl script to align RNA-seq single-end fastq files with [STAR](https://github.com/alexdobin/STAR) and then extract transcript counts with [featureCounts/subread](http://subread.sourceforge.net/).
This is designed to be a convenience for RNA-Seq processing, and to save time writing scripts.

## options

**-a** annotation file for featureCounts.  Accepts file afterward.  Optional.

**-d** debug.  Only print commands that would be executed.  Does not accept arguments.  Optional.

**-g** genome directory for STAR alignment.  Accepts directory as argument.  Required.

## examples

### specify an annotation file and a STAR genome directory

```shell
perl rnaseq.pl -a ~/GENE_DATA/gencode.v26.annotation.gtf -g ~/GENE_DATA/hg38
```
### only do debugging, just to see the commands that rnaseq.pl will output

```shell
perl rnaseq.pl -a  ~/GENE_DATA/gencode.v26.annotation.gtf -g ~/GENE_DATA/hg38 -d
```
