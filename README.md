# soma-cnv
Detect somatic copy number changes in low-depth sequencing data

## Overview

soma-cnv is a suite of tools to detect somatic copy number changes using low depth whole genome sequencing data.  Somatic copy number changes often occur in only a small proportion of the cells contributing to the sequenced DNA, and consequently manifest as subtle signals.  soma-cnv combines data from adjacent loci to increase its sensitivity to these subtle changes, trading positional resolution for sensitivity.  The output of soma-snv is a segmentation of the genome into copy number regions with estimated allele ploidies, and an estimation of the aneuploid fraction present in a sample.

![Example of soma-cnv detections](/docs/example.png?raw=true "Somatic copy number changes detected by soma-cnv")


## Workflow


## Preparation of whitelist and calibration files

soma-loh requires platform-specific calibration files for accurate CNV detection.  Whitelist and calibration files suitable for human data generated on the HiSeqX sequencer from TruSeq Nano libraries, and mapped to hs37d5 with PhiX decoy, are supplied in data/.  For a different species, platform, or pipeline, new files are likely to be required.  This section describes the creation of these files.

### 1. Locus whitelist
The first file that needs to be created is a set of whitelist loci to use for CNV determination.  These loci should be a set of positions on the genome which can be genotyped with very high reliability on the chosen platform.  The details of whitelist definition will differ between platforms, but in general whitelist loci should have the following properties:
* Consistently high genotyping rate across multiple samples.
* Consistent depth of sequencing across multiple samples.
* Autosomal, uniquely mapping, high mdust complexity, not in a repeat region.
* Biallelic
* High variant allele frequency (>= 5%) in the cohort, and no excess heterozygosity (eg judged by HWE test)

Once whitelist variants have been identified, they should be saved to a tab-delimited table with header, in the following format: ```chrom  pos  ref  alt```.  Example:
```
chrom	pos	ref	alt
1	100000012	G	T
1	100000827	C	T
1	100002882	T	G
1	100004726	G	A
1	100005477	G	A
1	100006117	G	A
1	100006734	C	T
```

### 2. GC content
After loci have been defined, GC content covariates need to be calculated at these loci for per-sample GC correction.  The script in util/generate_gc.py has been written for this purpose.  Run as follows:
```
python generate_gc.py <reference.fa> <whitelist.loci> > whitelist.gc.tsv
```

Where ```reference.fa``` is the reference fasta, ```whitelist.loci``` is the set of whitelist loci from step 1, and ```whitelist.gc``` is the gc covariate file.

### 3. Affinity correction



## Example

## Future