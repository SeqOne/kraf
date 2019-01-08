import docopt
import algorithm, strutils
import os
import hts
import tables

# KRAF stands for "K-mer Recalibrated Allele Frequency"

var
  fai : Fai
  fasta : string
  bam : string
  vcf : string
  t : bool
  k : int = 31
  padding : int = 5 # Padding around the deletion (sould be < k)

let doc = format("""
Usage: kraf <fasta> <bam> <vcf>

Arguments:

  <fasta>         Reference genome in FASTA format.
  <bam>           BAM files with aligned reads
  <vcf>           VCF files with called variants
""")

let args = docopt(doc)

fasta = $args["<fasta>"]
bam   = $args["<bam>"]
vcf   = $args["<vcf>"]

# Open fasta index to seek for k-mers
t = open(fai, fasta)

# Handy procedure to get sequence from 1-based input coordinates
proc getSeq(chr: string, start_pos: int, end_pos: int): string =
  get(fai, chr, start_pos - 1, end_pos - 1)

var
  comp_dna = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'a': 't',
    'c': 'g',
    'g': 'c',
    't': 'a',
    'N': 'N',
    'n': 'n',
  }.toTable

proc revcomp(s: var string) : string =
  result = newString(s.len)
  for i,c in s:
    result[s.high - i] = comp_dna[c]

proc canonicalKmer(s: var string) : string =
  var
    revcomp_s : string = revcomp(s)
  if s < revcomp_s:
    result = s
  else:
    result = revcomp_s

proc median(xs: seq[float]): float =
  var ys = xs
  sort(ys, system.cmp[float])
  0.5 * (ys[ys.high div 2] + ys[ys.len div 2]) 

proc kmersMedian(kmer_hash: var Table[string, int], kmers: var seq[string]) : float =
  var
    v: int
    values: seq[float] = @[]
  for kmer in kmers:
    v = kmer_hash[kmer]
    # Do not account for 0 count that could arise from an other variant beeing close to the deletion
    if v != 0:
      values.add(v.float)
  result = median(values)
  

# Set deletion position
var
  chr = "17"
  spos = 7573967 # This base is the first deleted base
  epos = 7574002 # This base is the last deleted base

# TODO Should we get the canonic k-mers for all k-mers that are going to be counted in the BAM ?

type
  KmerCount = tuple
    kmer: string
    count: int

var
  kleft : string
  kright : string
  kmer : string
  canon_kmer : string
  revcomp_kmer : string
  kmer_hash = initTable[string, int]() # Hash table to count for k-mers
  del_kmers : seq[string] = @[]
  ref_kmers : seq[string] = @[]

# Build k-mers for the deletions
for i in (spos - k + padding)..(spos - padding):
  # Get left part of the k-mer
  kleft = getSeq(chr, i, spos - 1)
  # Get right part of the k-mer
  kright = getSeq(chr, epos + 1, epos + (k - len(kleft)))
  # Merge both part of the k-mer targeting the deletion
  kmer = kleft & kright
  # Get canonical k-mer
  kmer = canonicalKmer(kmer)
  del_kmers.add(kmer)
  # Init the k-mer hash for this k-mer
  kmer_hash[kmer] = 0

# Build k-mers for the reference (not deleted) sequence
for i in (spos - k + padding)..(epos - padding + 1):
  kmer = getSeq(chr, i, i + k - 1)
  # Get canonical k-mer
  kmer = canonicalKmer(kmer)
  ref_kmers.add(kmer)
  # Init the k-mer hash for this k-mer
  kmer_hash[kmer] = 0

# stderr.writeLine '------------\nDeletion k-mers\n'
# for k in del_kmers:
#   stderr.writeLine k

# stderr.writeLine '\n------------\nReference k-mers\n'
# for k in ref_kmers:
#   stderr.writeLine k

# TODO Before counting k-mers in the BAM we should check the reference 
# to remove k-mers with duplicated location on the genome, etc.

# Loop over the BAM files to count the k-mers
# open a bam and look for the index.
var 
  b:Bam
  read: string
  revcomp_read: string
  nb_reads: int = 0

# Open the bam file
open(b, bam, index=true)

echo "Reading read sequence in BAM file to recalibrate VAF"
for record in b:
  # Skip duplicate reads based on flag
  if record.flag.dup:
    continue
  # TODO skip bad sequence based on quality
  # TODO use canonical k-mers instead of counting both
  read = sequence(record, read)
  revcomp_read = revcomp(read)
  nb_reads += 1
  if (nb_reads mod 10000 == 0):
    echo nb_reads, " parsed ..."

  for i in 0..(len(read) - k):
    kmer = read[i .. i + k - 1]
    revcomp_kmer = revcomp_read[i .. i + k - 1]
    if kmer < revcomp_kmer:
      canon_kmer = kmer
    else:
      canon_kmer = revcomp_kmer
    if kmer_hash.hasKey(canon_kmer):
      kmer_hash[canon_kmer] += 1

  # # Count for Revcomp k-mers
  # read = revcomp(read)
  # for i in 0..(len(read) - k):
  #   kmer = read[i .. i + k - 1]
  #   if kmer_hash.hasKey(kmer):
  #     kmer_hash[kmer] += 1

var
  ref_median: float = kmersMedian(kmer_hash, ref_kmers)
  alt_median: float = kmersMedian(kmer_hash, del_kmers)
  refined_vaf = alt_median / ref_median

echo "Ref median: ", ref_median
echo "Alt median: ", alt_median
echo "Refined VAF: ", refined_vaf

# TODO Only replace original VAF at certain condition:
#         - New VAF is higher than the previous one
#         - New VAF is no more than X fold higer
# TODO Keep the old VAF in a new field : OLD_AF