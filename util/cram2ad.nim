let doc = """
Convert CRAM mapped reads to an allele depth file for use by soma-cnv.

Usage:
  cram2ad [options] <loci.tsv> <reference.fa> <reads.cram>
  cram2ad (-h | --help)
  cram2ad --version

Options:
  -h --help      Show this screen.
  --version      Show version.
  --minbq=<n>    Minimum base quality [default: 30].
  --minmq=<n>    Minimum mapping quality [default: 40].
  --threads=<n>  Number of CRAM decoding threads [default: 0].

The input CRAM file should have duplicates marked.

Mark Pinese (m.pinese@garvan.org.au)
12 April 2019
"""

import docopt
import hts
import os
import strutils

const CRAM_ALL_FIELDS = int(8191)
const CRAM_DISCARD_FIELDS = SAM_QNAME.int or SAM_RNAME.int or SAM_AUX.int or SAM_RGAUX.int or SAM_RNEXT.int or SAM_PNEXT.int or SAM_TLEN.int
const CRAM_WANTED_FIELDS = CRAM_ALL_FIELDS - CRAM_DISCARD_FIELDS


proc quickpile(bam: Bam, chrom: string, position: int, ref_allele: char, alt_allele: char, min_mapq: uint8, min_baseq: uint8): tuple[nref: int, nalt: int, noth: int] =
  # Adapted from CIGAR skipping code by Brent Pedersen (https://brentp.github.io/post/no-pile/)
  # position is 1-based inclusive
  var
    nalt = 0
    nref = 0
    noth = 0

  for aln in bam.query(chrom, position - 1, position):
    var
      off = aln.start
      qoff = 0
      roff_only = 0

    # Skip reads with unwanted flags
    if aln.flag.unmapped or aln.flag.secondary or aln.flag.qcfail or aln.flag.dup:
      continue

    # Skip reads with too low MQ
    if aln.mapping_quality < min_mapq:
      continue

    for event in aln.cigar:
      let cons = event.consumes
      if cons.query:
        qoff += event.len
      if cons.reference:
        off += event.len
        if not cons.query:
          roff_only += event.len

      # stop once we've parsed enough CIGAR to infer the
      # target locus's position in the query (read)
      if off > position:
        break

    # since each cigar op can consume many bases
    # calc how far past the requested position
    let
      over = off - position - roff_only
      qpos = qoff - over - 1

    if qpos == -1:
      continue

    # verify the quality at this base
    let qual = aln.base_quality_at(qpos)
    if qual >= min_baseq:
      # get the base 
      let base = aln.base_at(qpos)
      if base == ref_allele:
        nref += 1
      elif base == alt_allele:
        nalt += 1
      else:
        noth += 1

  return (nref, nalt, noth)


proc main() =
  let
    args = docopt(doc, version = "cram2ad 1.0.0")
    variants_file = open($args["<loci.tsv>"])
    reference_path = $args["<reference.fa>"]
    cram_path = $args["<reads.cram>"]
    min_mapq = ($args["--minmq"]).parseInt.uint8
    min_baseq = ($args["--minbq"]).parseInt.uint8
    threads = ($args["--threads"]).parseInt

  var
    cram: Bam
    firstline = true

  open(cram, cram_path, index=true, fai=reference_path, threads=threads)
  if cram.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, CRAM_WANTED_FIELDS) != 0:
    quit "Error setting CRAM decode options"

  for line in variants_file.lines:
    if firstline:
      firstline = false
      continue
    var
      lineparts = line.split('\t')
      nref: int
      nalt: int
      noth: int
    let
      chrom = lineparts[0]
      pos = parseInt(lineparts[1])
      ref_allele = lineparts[2][0]
      alt_allele = lineparts[3][0]
    (nref, nalt, noth) = quickpile(cram, chrom, pos, ref_allele, alt_allele, min_mapq, min_baseq)
    if nref == 0 or nalt == 0:
      continue
    echo chrom, '\t', pos, '\t', nref + nalt, '\t', nalt


main()
