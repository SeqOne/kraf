import docopt
import algorithm, strutils
import os
import hts
import tables

# KRAF stands for "K-mer Recalibrated Allele Frequency"

const kraf_version = "0.1.0"

type
  KmerSeq = seq[string]
  KmerArray = seq[KmerSeq]
  VarKmers = tuple
    ref_kmers: KmerArray
    alt_kmers: KmerArray
    var_key: string

# Handy procedure to get sequence from 1-based input coordinates
proc getSeq(fai: Fai, chr: string, start_pos: int, end_pos: int): string =
  get(fai, chr, start_pos - 1, end_pos - 1).toUpperAscii()

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

# TODO Create procedures to generate the k-mers for SNV, Indel, Ins & Del event
# This will break down the code into smaller pieces mores testable and will
# prepare the ground for Kraf evolutions.
# proc deletionKmers(chr: var string, spos: int, epos: int, k: int, padding: int) : KmerSeq =

proc median(xs: seq[float]): float =
  var ys = xs
  sort(ys, system.cmp[float])
  0.5 * (ys[ys.high div 2] + ys[ys.len div 2])

proc kmersMedian(kmer_hash: var Table[string, int], kmers: var KmerArray) : float =
  var
    v: int
    values: seq[float] = @[]
  for kmer_seq in kmers:
    v = 0
    # For all k-mers at this position we merge their values
    for kmer in kmer_seq:
      if kmer_hash.hasKey(kmer):
        v += kmer_hash[kmer]
        # Do not account for 0 count that could arise from an other variant beeing close to the deletion
      else:
        stderr.writeLine "Kmer not in hash table !!!"
    if v != 0:
      values.add(v.float)
  stderr.writeLine repr(kmers)
  stderr.writeLine repr(values)
  if len(values) == 0:
    result = 0
  else:
    stderr.writeLine repr(values)
    result = median(values)

var
  fai : Fai
  v : VCF
  v_query : VCF
  b : Bam
  fasta : string
  bam : string
  vcf : string
  t : bool
  verbose : bool = false
  k : int = 31
  padding : int = 10 # Padding around the deletion (sould be < k)
  threads : int = 4
  min_del_size : int = 10

let doc = format("""
version: $version

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
""", @["version", kraf_version])

let args = docopt(doc)

# Parse arguments
fasta         = $args["<fasta>"]
bam           = $args["<bam>"]
vcf           = $args["<vcf>"]

# Parse options
k             = parseInt($args["-k"])
padding       = parseInt($args["-p"])
threads       = parseInt($args["-t"])
min_del_size  = parseInt($args["--min-del"])

if args["--verbose"]:
  verbose = true

if padding >= k:
  quit "Padding value must be inferior to k"

# Open fasta index to seek for k-mers
if not open(fai, fasta):
  quit "Failed to open FASTA file"

var
  chr : string
  spos : int # This base is the first deleted base
  epos : int # This base is the last deleted base
  kleft : string
  kright : string
  kmer : string
  var_key : string
  var_key2 : string
  canon_kmer : string
  revcomp_kmer : string
  kmer_hash = initTable[string, int]() # Hash table to count for k-mers
  variant_hash = initTable[string, VarKmers]()
  vark : VarKmers

# Parse VCF to construct k-mers for recalibrating VAF of deletions
# Open VCF file
if not open(v, vcf):
  quit "Failed to open VCF file"

# Open a second handler on the VCF to query on intervals for bi-allelic variants
if not open(v_query, vcf):
  quit "Failed to open VCF file"

# Quit if VCF has more than one sample!!!
if v.n_samples != 1:
  quit "only works for 1 sample"

for rec in v:
  # Test if we have a deletion (strict mode for now, no indel style)
  # TODO Adapt Kraf for delins and insertions
  if len(rec.REF) > min_del_size:
    for alt in rec.ALT:
      if len(alt) == 1:
        chr = $rec.CHROM
        spos = rec.POS + 1
        epos = rec.POS + len(rec.REF) - 1

        var_key = chr & '-' & $rec.POS & '-' & rec.REF & '-' & alt

        if verbose:
          stderr.writeLine "Adding deletion : ", var_key

        # Init the Variant Tuple
        vark = (ref_kmers: @[], alt_kmers: @[], var_key: var_key)

        # Build k-mers for the deletions
        for i in (spos - k + padding)..(spos - padding):
          # Get left part of the k-mer
          kleft = getSeq(fai, chr, i, spos - 1)
          # Get right part of the k-mer
          kright = getSeq(fai, chr, epos + 1, epos + (k - len(kleft)))
          # Merge both part of the k-mer targeting the deletion
          kmer = kleft & kright
          # Get canonical k-mer
          kmer = canonicalKmer(kmer)
          # Add k-mer to the sequence of k-mers at this position
          var a : seq[string ]= @[]
          a.add(kmer)
          vark.alt_kmers.add(a)
          # Init the k-mer hash for this k-mer
          kmer_hash[kmer] = 0

        # Build k-mers for the reference (not deleted) sequence
        var
          ref_start = spos - k + padding # First k-mer position
          ref_end = epos - padding + 1 # Last k-mer position
          alt_window_end = ref_end + k - 1 # Last base of the last k-mer

        
        for i in ref_start..ref_end:
          kmer = getSeq(fai, chr, i, i + k - 1)
          # Get canonical k-mer
          kmer = canonicalKmer(kmer)
          var a : seq[string ]= @[]
          a.add(kmer)
          vark.ref_kmers.add(a) 
          # Init the k-mer hash for this k-mer
          kmer_hash[kmer] = 0

        # Build k-mers for bi-allelic variants
        for rec in v_query.query(chr & ':' & $ref_start & '-' & $alt_window_end):
          var ref_len = len(rec.REF)
          for alt in rec.ALT:

            # If this is the current variant, we skip it
            var_key2 = chr & '-' & $rec.POS & '-' & rec.REF & '-' & alt

            # We do not want to add the current variant as it's own co-occurence
            if var_key == var_key2:
              continue

            if verbose:
              stderr.writeLine "  - Adding co-occurence variant" & var_key2
            
            var 
              alt_len = len(alt)
              ins_len = alt_len - 1
              ins_seq = alt[1 .. ins_len]
              ins_pos = rec.POS + 1
              del_len = ref_len - 1
              del_pos = rec.POS + 1
              del_end = del_pos + del_len - 1

            # SNV case
            if ref_len == alt_len and ref_len == 1:
              #echo rec
              for i in max(ref_start,rec.POS - k + 1)..min(ref_end,rec.POS):
                # Get the position of the SNP in the k-mer
                var snp_kmer_pos = rec.POS - i
                # Get reference sequence at this position
                # kmer = getSeq(fai, chr, i, i + k - 1)
                kmer = vark.ref_kmers[i - ref_start][0]
                #echo kmer
                # Update the kmer with alternate base
                kmer[snp_kmer_pos] = alt[0]
                kmer = canonicalKmer(kmer)
                #echo kmer
                # Add the kmer to the kmer array at this position
                vark.ref_kmers[i - ref_start].add(kmer)
                kmer_hash[kmer] = 0
            # Insertion case
            elif ref_len == 1 and alt_len > 1:
              # If the insertion can be targeted by a k-mer
              if ins_len + 2 * padding <= k:
                for i in max(ref_start,ins_pos - k + padding + ins_len)..min(ref_end,ins_pos - padding):
                  # Get left part of the k-mer
                  kleft = getSeq(fai, chr, i, rec.POS)
                  # Get right part of the k-mer
                  kright = getSeq(fai, chr, rec.POS + 1, rec.POS + (k - ins_len - len(kleft)))
                  # Merge both part of the k-mer targeting the insertion 
                  kmer = kleft & ins_seq & kright
                  # Get canonical k-mer
                  kmer = canonicalKmer(kmer)
                  # Add the kmer to the kmer array at this position
                  vark.ref_kmers[i - ref_start].add(kmer)
                  # Init the k-mer hash for this k-mer
                  kmer_hash[kmer] = 0
              else:
                # TODO In that case, maybe we should not recalibrate de VAF as we will
                # overstimate its value ...
                stderr.writeLine "Skip insertion variant is too long for k-mer tagging"
            # Deletion case
            elif alt_len == 1 and ref_len > 1:
              for i in max(ref_start, del_pos - k + padding)..min(ref_end, del_pos - padding):
                # Get left part of the k-mer
                kleft = getSeq(fai, chr, i, rec.POS)
                # Get right part of the k-mer
                kright = getSeq(fai, chr, del_end + 1, del_end + (k - len(kleft)))
                # Merge both part of the k-mer targeting the deletion
                kmer = kleft & kright
                # Get canonical k-mer
                kmer = canonicalKmer(kmer)
                # Add the kmer to the kmer array at this position
                vark.ref_kmers[i - ref_start].add(kmer)
                # Init the k-mer hash for this k-mer
                kmer_hash[kmer] = 0
            # Case for indel
            else:
              # If the insertion can be targeted by a k-mer
              if ins_len + 2 * padding <= k:
                for i in max(ref_start,del_pos - k + padding + ins_len)..min(ref_end,del_pos - padding):
                  # Get left part of the k-mer
                  kleft = getSeq(fai, chr, i, rec.POS)
                  # Get right part of the k-mer
                  kright = getSeq(fai, chr, del_end + 1, del_end + (k - ins_len - len(kleft)))
                  # Merge both part of the k-mer targeting the insertion 
                  kmer = kleft & ins_seq & kright
                  # Get canonical k-mer
                  kmer = canonicalKmer(kmer)
                  # Add the kmer to the kmer array at this position
                  vark.ref_kmers[i - ref_start].add(kmer)
                  # Init the k-mer hash for this k-mer
                  kmer_hash[kmer] = 0
              else:
                stderr.writeLine "Skip delins variant is too long for k-mer tagging"

        # Add the variant to the variant hash
        variant_hash[var_key] = vark

close(v)

# TODO Before counting k-mers in the BAM we should check the reference
# to remove k-mers with duplicated location on the genome, etc.

# Loop over the BAM files to count the k-mers
# open a bam and look for the index.
var
  read: string
  revcomp_read: string
  nb_reads: int = 0

# Open the bam file
open(b, bam, index=true, threads=threads)

if b == nil:
    quit "could not open bam file"
if b.idx == nil:
    quit "could not open bam index"

if verbose:
  stderr.writeLine "Reading read sequence in BAM file to recalibrate VAF"

for record in b:
  # Skip duplicate reads based on flag
  if record.flag.dup:
    continue
  # TODO skip bad sequence based on quality
  read = sequence(record, read)
  #revcomp_read = revcomp(read)
  nb_reads += 1
  if verbose and nb_reads mod 10000 == 0:
    stderr.writeLine nb_reads, " parsed ..."

  for i in 0..(len(read) - k):
    kmer = read[i .. i + k - 1]
    kmer = canonicalKmer(kmer)
    # revcomp_kmer = revcomp_read[i .. i + k - 1]
    # if kmer < revcomp_kmer:
    #   canon_kmer = kmer
    # else:
    #   canon_kmer = revcomp_kmer
    if kmer_hash.hasKey(kmer):
      kmer_hash[kmer] += 1

# Now we open the VCF a second time to print updated VAF
var
  ref_median : float
  alt_median : float
  depth : int32
  alternate_obs : int32
  reference_obs : int32
  refined_vaf : float32
  old_vaf : float32

t = open(v, vcf)

# Add headers for AO and AD (if not yet present)
# ##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
if v.header.add_format("AO", "A", "Integer", "Count of full observations of this alternate haplotype.") != Status.OK:
  quit "unable to add AO to the header"
# ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
if v.header.add_format("AD", "R", "Integer", "Allelic depths for the ref and alt alleles in the order listed.") != Status.OK:
  quit "unable to add AD to the header"
# ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
if v.header.add_format("AF", "1", "Float", "Allele Frequency.") != Status.OK:
  quit "unable to add AF to the header"

if v.header.add_format("OLD_AF", "1", "Float", "Old AF value that have been replaced by KRAF") != Status.OK:
  quit "unable to add OLD_AF to the header"
if v.header.add_format("OLD_AO", "1", "Float", "Old AO value that have been replaced by KRAF") != Status.OK:
  quit "unable to add OLD_AF to the header"
if v.header.add_format("OLD_AD", "R", "Integer", "Old AD value that have been replaced by KRAF") != Status.OK:
  quit "unable to add OLD_AF to the header"

# Print headers
stdout.write v.header

for rec in v:
  for alt in rec.ALT:
    var_key = $rec.CHROM & '-' & $rec.POS & '-' & rec.REF & '-' & alt
    if variant_hash.hasKey(var_key) :

      # Get the deletion ref & alt k-mers and compute median counts for both
      vark = variant_hash[var_key]
      ref_median  = kmersMedian(kmer_hash, vark.ref_kmers)
      alt_median  = kmersMedian(kmer_hash, vark.alt_kmers)

      # We failed at finding k-mers for either the reference or the alternative
      # We therefor skip the recalibration
      if (ref_median + alt_median) == 0:
        continue

      # Compute the refined VAF
      refined_vaf = alt_median / (ref_median + alt_median)

      var floats = newSeq[float32](1)
      var one_int = newSeq[int32](1)
      var two_ints = newSeq[int32](2)
      var two_ints_ad_field = newSeq[int32](2)

      if rec.format.get("DP", one_int) == Status.OK:
        depth = one_int[0]
      else:
        quit "missing DP field for a deletion variant"

      if rec.format.get("AD", two_ints_ad_field) == Status.OK:
        old_vaf = two_ints_ad_field[1] / depth
      elif rec.format.get("AF", floats) == Status.OK:
        old_vaf = floats[0]
      else:
        quit "missing AF or AD field to compute previous VAF"

      if verbose:
        stderr.writeLine "\nUpdated DELETION: ", var_key
        stderr.writeLine "Ref median: ", ref_median
        stderr.writeLine "Alt median: ", alt_median
        stderr.writeLine "Refined VAF: ", refined_vaf
        stderr.writeLine "Old VAF: ", old_vaf

      # Skip refined VAF if is lower than the previous one
      if refined_vaf <= old_vaf:
        continue
      elif refined_vaf > 1:
        continue
      # elif refined_vaf > 4 * old_vaf: # This seems very unlikely to happened and be true
      #   continue

      # Keep old values for AF, AO & AD if they are found
      if rec.format.get("AF", floats) == Status.OK:
        if rec.format.set("OLD_AF", floats) != Status.OK:
          quit "error setting OLD_AF in VCF"

      if rec.format.get("AO", one_int) == Status.OK:
        if rec.format.set("OLD_AO", one_int) != Status.OK:
          quit "error setting OLD_AO in VCF"

      if rec.format.get("AD", two_ints) == Status.OK:
        if rec.format.set("OLD_AD", two_ints) != Status.OK:
          quit "error setting OLD_AD in VCF"

      # We compute Alternate observation from original DP as we could
      # have underestimated depth with exact k-mers due to sequencing errors
      alternate_obs = (refined_vaf * depth.float()).int32

      # Set AF
      floats[0] = refined_vaf
      if rec.format.set("AF", floats) != Status.OK:
        quit "error setting AF in VCF"

      # Set AO
      one_int[0] = alternate_obs
      if rec.format.set("AO", one_int) != Status.OK:
        quit "error setting AO in VCF"

      # Set AD
      two_ints[0] = depth - alternate_obs
      two_ints[1] = alternate_obs
      if rec.format.set("AD", two_ints) != Status.OK:
        quit "error setting AD in VCF"

  # Echo the record (modified or not)
  stdout.write rec.tostring()