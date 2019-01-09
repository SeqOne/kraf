# Kraf (Pronunce Carafe) - Kmer based Recalibration of Allele Frequency

Allele frequency estimation for long deletions made by variant caller (eg. freebayes) are usually underestimated
due to the fact that these softwares only compare reads with identified breakpoint to
the coverage at the position of the deletion. Also some read may also contained the deletion,
but near their extremity, and the alignement of the deletion was not possible resulting
in softclipping.

## Method

Kraf uses a k-mers strategy to recalibrate the Allele frequency of these deletions
using the signal observed on the k-mer level. This strategy is alignment-free and
results for this variants into a better estimation of the AF.

## Usage

Kraf takes in input the reference genome in fasta format, a BAM file and the VCF with
called variants. It output on STDOUT a modified VCF with updated VAF for deletions.

```
Usage: kraf [options] <fasta> <bam> <vcf>

Arguments:
    <fasta>           Reference genome in FASTA format.
    <bam>             BAM files with aligned reads
    <vcf>             VCF files with called variants

Options:
    -k <int>          K-mer size to use for recalibration. [default: 31]
    -p <int>          Padding size of k-mers around the deletion. [default: 5]
    -t <int>          Number of decompression threads. [default: 4]
    --min-del <int>   Minimun deletion size elligible for AF recalibration [default: 10]
    -v, --verbose     Enable verbose mode
```

KRAF will set the following values in the FORMAT field of the sample :
- AO : Number of alternate observations
- AD : Allelic depths for the ref and alt alleles in the order listed.
- AF : Allele Frequency

When an AF, AO and AD are beeing replaced in the FORMAT field a new TAG named "OLD_TAG" (ex: OLD_AF) 
is appended in order to keep the trace of the previous values computed by the caller and replaced
by KRAF.

Beware of these Kraf limitations :

- Only works for 1 sample VCF file
- Only deals with deletions

## Install

First install dependencies : nim compilator and hts-nim library

```
apt-get install nim
nimble install -y hts
```

Compile kraf binary and place it somewhere available from your `$PATH`:

```
git clone https://gitlab.seq.one/workset/kraf.git && cd kraf
make
```
