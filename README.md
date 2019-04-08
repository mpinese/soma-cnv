# soma-cnv
Detect somatic copy number changes in low-depth sequencing data

## Overview

soma-cnv is a suite of tools to detect somatic copy number changes using low depth whole genome sequencing data.  Somatic copy number changes often occur in only a small proportion of the cells contributing to the sequenced DNA, and consequently manifest as subtle signals.  soma-cnv combines data from adjacent loci to increase its sensitivity to these subtle changes, trading positional resolution for sensitivity.  The output of soma-snv is a segmentation of the genome into copy number regions with estimated allele ploidies, and an estimation of the aneuploid fraction present in a sample.

Examples of somatic copy number changes detected by soma-cnv:

![Example of soma-cnv detections](/docs/examples_1_2.png?raw=true "Somatic copy number changes detected by soma-cnv")

soma-cnv's key properties:
* Highly sensitive: calling copy number variation in low purity samples (down to ~ 5% aneuploid nucleated cell fraction)
* Suitable for low-depth data: soma-cnv was developed to work on human WGS intended for germline genotyping (30-40X mean depth).
* Approximate resolution of CNV only: soma-cnv's spatial resolution depends on the loci it is supplied.  Typically CNV is localised to within 10 kb, but soma-cnv does not give base-level breakpoint resolution.
* Assumes the presence of a single clonal aneuploid fraction mixed among diploid cells.  Deviations from this assumption (eg aneuploid subclones) will degrade performance.
* Works on a per-sample basis.
* Requires calibration files.  One set of calibration files is supplied, but it may be necessary to generate calibration files for different experimental setups.  Although soma-cnv's primary calling workflow works in single-sample mode, the generation of calibration files requires a large collection of platform-matched samples.


## Requirements

* Mapped sequencing data
* A variant caller (eg GATK HaplotypeCaller, samtools.  Tested with GATK 3.7-0-gcfedb67)
* R (tested with v3.4.2)
  * mgcv (tested with 1.8-22)
  * plyr (tested with 1.8.4)
  * ggplot2 (tested with 2.2.1)
  * tvd (tested with 0.1.0)


## Workflow

The main soma-snv workflow (steps 2 and 3) operates on a per-sample basis.  If the generation of platform calibration files is required (step 1), this requires the analysis of many platform-matched samples in parallel.

1. (optional) Prepare platform-specific calibration files.
2. For each sample, collect allele-specific depths at whitelist loci.
3. For each sample, for a soma-cnv model and identify aneuploid regions.

## Example

### 1. Cohort: prepare whitelist and calibration files

soma-cnv requires platform-specific calibration files for accurate CNV detection.  Whitelist and calibration files suitable for human data generated on the HiSeqX sequencer from TruSeq Nano libraries, and mapped to hs37d5 with PhiX decoy, are supplied in data/.  For a different species, platform, or pipeline, new files are likely to be required.  This section describes the creation of these files.

Although soma-cnv runs in a single-sample mode, the generation of calibration files requires sequencing data from many (ideally over 100) platform-matched samples.  These samples must be from different individuals; the use of technical replicates to generate calibration data will reduce sensitivity.


#### A. Locus whitelist
The first file that needs to be created is a set of whitelist loci to use for CNV determination.  These loci should be a set of positions on the genome which can be genotyped with very high reliability on the chosen platform.  The details of whitelist definition will differ between platforms, but in general whitelist loci should have the following properties:
* Consistently high genotyping rate across multiple samples.
* Consistent and typical depth of sequencing across multiple samples.
* Autosomal, uniquely mapping, high mdust complexity, not in a repeat region.
* Intermediate GC content, for example in [0.3, 0.55] in a 100 bp window around the locus.
* Biallelic
* High variant allele frequency (>= 5%) in the cohort, and no excess heterozygosity (eg judged by HWE test)

Once whitelist variants have been identified, they should be saved to a tab-delimited table with header, in the following format: `chrom	pos	ref	alt`.  Example:
```
chrom	pos	ref	alt
1	100000012	G	T
1	100000827	C	T
1	100002882	T	G
1	100004726	G	A
1	100005477	G	A
```

#### B. GC content
After loci have been defined, GC content covariates need to be calculated at these loci for per-sample GC correction.  The script in util/generate_gc.py has been written for this purpose.  Run as follows:
```
python generate_gc.py reference.fa whitelist.loci.tsv > whitelist.gc.tsv
```

Where `reference.fa` is the reference fasta, `whitelist.loci.tsv` is the set of whitelist loci from step A, and `whitelist.gc.tsv` is the generated gc covariate file.


#### C. Affinity calculation
The final calibration file required is the locus affinity file `whitelist.affinity.tsv`, which is used to correct for platform-specific sequencing depth differences between loci.  The affinity file is a tab-separated table with header `chrom	pos	affinity`, for example:
```
chrom	pos	affinity
1	100000012	0.971763973729247
1	100000827	1.03221998575009
1	100002882	1.0519381589067
1	100004726	0.981624465488513
1	100005477	1.03637637829641
```
`chrom` and `pos` correspond to the loci in `whitelist.loci.tsv`, and must contain the same loci in the same order.  `affinity` is the mean normalised depth at `chrom:pos` across the calibration cohort.  To compute `affinity` values:
1. Collect a cohort of calibration samples.  These samples should all be technically similar: eg matched for DNA source, extraction method, sequencing, and mapping, and be the same as the test samples for these metrics.  With a sufficiently large cohort it is acceptable to use the test samples as calibration samples also.
2. For each calibration sample i of n total:
   1. At each locus j of m total loci in `whitelist.loci.tsv`, measure the duplicate-excluded sequencing depth ![d_i_j](/docs/eqn_d_i_j.gif?raw=true "d_i_j").
   2. Calculate the mean depth across all loci in `whitelist.loci.tsv`, ![eqn_d_bar_i](/docs/eqn_d_bar_i.gif?raw=true "eqn_d_bar_i").
   3. At each locus in `whitelist.loci.tsv`, calculate the normalised depth, ![eqn_c_i_j](/docs/eqn_c_i_j.gif?raw=true "eqn_c_i_j").
3. For each locus, calculate the mean normalised depth across all samples, ![eqn_c_j](/docs/eqn_c_j.gif?raw=true "eqn_c_j").
4. Normalise the mean normalised depth to have mean of 1, ![eqn_a_hat_j](/docs/eqn_a_hat_j.gif?raw=true "eqn_a_hat_j").

The resultant ![a_hat_j](/docs/eqn_hat_a_j.gif?raw=true "\hat{a_j}") are the affinity values to insert into `whitelist.affinity.tsv`.

Note that if your data are particularly variable, a robust alternative (eg median) to the means in steps 2i and 3 will also work.  The normalisation in step 4 should remain the arithmetic mean.


### 2. Per-sample: collection of allele-specific depths

This step examines the input mapped BAM file at the loci given in `whitelist.loci.tsv`, and reports the total and alt allele depths at each heterozygous locus in a tab-delimited table with format `<chrom>	<pos>	<dp>	<ad>`
Note that this file has no header, and that homozygous reference or alternate loci are not included.  Example:
```
1	943687	48	18
1	944564	38	22
1	1031540	33	22
1	1033999	39	22
1	1065296	31	20
```

This file can be generated by a number of methods, and soma-cnv should be quite robust to the choice.  Example code for GATK HaplotypeCaller:
```bash
# Convert the tsv-format whitelist loci into an interval_list for use by GATK.
xz -dc data/truseq_nano_hiseqX_hs37d5x.loci.tsv.xz | \
  awk '{print $1 ":" $2 "-" $2}' \
  > truseq_nano_hiseqX_hs37d5x.loci.interval_list

# Run GATK HC
# -ip 100 instructs HC to consider a region of 100 bp around each locus, to enable
# local haplotype reassembly.  Note that because of this, some additional variant
# loci may be reported (not just those in the interval_list), but these will be
# removed at the later R stage.
java -Xmx2G -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R <reference.fa> \
  -L truseq_nano_hiseqX_hs37d5x.loci.interval_list -ip 100 \
  -I sampleID.bam -o sampleID.temp.vcf

# Post-process the GATK VCF: extract het SNP loci and report depths.
awk -f util/filter_hc_vcf.awk < sampleID.temp.vcf \
  | xz -c \
  > sampleID.soma-cnv.hetdp.xz
```

Example for bcftools:
```bash
xz -dc data/truseq_nano_hiseqX_hs37d5x.loci.tsv.xz | \
  awk 'BEGIN {OFS="\t"} {print $1, $2-1, $2}' \
  > truseq_nano_hiseqX_hs37d5x.loci.bed

bcftools mpileup -f <reference.fa> -B -C 50 -q 40 -Q 30 -I -R truseq_nano_hiseqX_hs37d5x.loci.bed -a 'FORMAT/AD' sampleID.bam | \
  python3 utils/bcftools_vcf2ad.py data/truseq_nano_hiseqX_hs37d5x.loci.tsv.xz | \
  xz -c \
  > sampleID.soma-cnv.hetdp.xz
```


### 3. Per-sample: Fit soma-cnv model and identify aneuploid regions

```bash
Rscript soma-cnv.R \
  data/truseq_nano_hiseqX_hs37d5x.affinity.tsv \
  data/truseq_nano_hiseqX_hs37d5x.gc.tsv \
  sampleID.soma-cnv.hetdp.xz sampleID.soma-cnv.rds
```

The primary output is a RDS containing a single list.  This list has members:
* data: Data frame of input allele depth data, with additional annotation fields added by soma-cnv.
* opts: List of options supplied to the algorithm.
* fit: List with elements describing the results of the soma-cnv fit. Has members:
  * models: Data frame summarising the models tested and their fits to the data.
  * model_search: Data frame giving the parameter space searched to find the best model.
  * fit.orig: The best-fitting segmentation before segment merging.
  * fit: The best-fitting segmentation after segment merging.

fit$fit is the primary output of a soma-cnv run.  Each row of fit$fit describes a segment of the genome of consistent copy number.  The most relevant fields of fit$fit are:
* chrom, start_pos, end_pos: genomic coordinates of the segment.  1-based, inclusive.
* fit.k1, fit.k2: copy numbers for the two chromatids in this segment within the aneuploid cells.
* fit.f: the estimated cellular fraction of the sample that is aneuploid.  Note that this is always identical across rows.

For example, this excerpt of a fit$fit data frame (from the example below, unnecessary columns dropped) describes a complex event on chromosome 13 affecting RB1 (13:48303775-48481986), that is present in 36% of nucleated cells contributing to the sample.
```
chrom start_pos   end_pos fit.k1 fit.k2     fit.f
   13  19509358  40757910      1      1 0.3619173   <-- diploid (ie no CNV)
   13  40758208  42306139      1      0 0.3619173     <-- single copy loss of 1.5 Mb
   13  42307184  42708231      1      1 0.3619173   <-- diploid (ie no CNV)
   13  42711664  42886868      1      0 0.3619173     <-- single copy loss of 0.2 Mb
   13  42890109  46649721      1      1 0.3619173   <-- diploid (ie no CNV)
   13  46649723  50646218      1      0 0.3619173     <-- single copy loss of 4.0 Mb, including all of RB1.
   13  50656582  51330918      0      0 0.3619173       <-- biallelic loss of 0.7 Mb
   13  51331121  51872361      1      0 0.3619173     <-- single copy loss of 0.5 Mb
   13  51880895 115044332      1      1 0.3619173   <-- diploid (ie no CNV)
```

**Note that soma-cnv does not phase segments** -- for example, in the above, it is not implied that the single-copy deletion events on lines 2 and 4 affect the same molecule.


#### Diagnostics

soma-cnv optionally emits a PDF of diagnostic plots, which can be useful to identify fit issues.

##### Diagnostic plot 1: homozygote exclusion

![Diagnostic plot 1: homozygote exclusion](/docs/example2-diagnostics-pg1.jpg?raw=true "Diagnostic plot 1: homozygote exclusion")

This plot is used to verify that soma-cnv is only modelling data from heterozygous loci.  Shown above is an example of a good plot, in which the majority of data points are from het loci (blue), and only a few are from homozygous loci (red).  If a large number of data points have a VAF close to zero or one, and especially if these points have been called as heterozygous by soma-cnv, the model will likely not fit well.  Reexamination of the variant calling pipeline to ensure that homozygous loci aren't reported, or improved filtering of the whitelist loci, might help.


##### Diagnostic plot 2: data points and fit

![Diagnostic plot 2: data points and fit](/docs/example2-diagnostics-pg3.jpg?raw=true "Diagnostic plot 2: data points and fit")

This plot shows the overall data distribution and fit across the genome.  This example shows a good result, in which the majority of the genome is fit as diploid (green line), and a simple smooth across the genome agrees with this assignation (red line in upper plot, largely obscured by green line).  In addition, the VAF plot is consistently clustered around 0.5, with no indication of contamination.  Failure modes are revealed on this plot by the following:
* Highly noisy depth data (red line in upper plot highly variable, often accompanied by a large number of aneuploid segment calls)
* Contamination (a stripe of points is visible near VAF = 0 and VAF = 1 on the lower plot, often accompanied by the green line in the lower plot being split across the whole genome).
* Overdispersion (the green line in the lower plot is split across the whole genome, without a stripe of points visible near VAF = 0 and VAF = 1).

Noisy depth data and overdispersion may indicate an issue with calibration, especially if this is observed in the majority of samples.  In that case, generation of custom calibration data may resolve the problem.  Contamination cannot be addressed except by resequencing, ideally from a freshly-collected sample.


#### Example exploration of the output RDS.

```R
result = readRDS("docs/example.rds")

names(result)
# [1] "data" "fit"  "opts"

head(result$data)
#        chrom     pos dp ad  affinity gc100 gc200  gc400     gc600   gc800 pois.lambda         prRR         prAA  het
# 289282     1  943687 48 18 0.9622546  0.48 0.445 0.4275 0.4633333 0.48375    40.33837 3.470682e-20 1.321768e-55 TRUE
# 289493     1  944564 38 22 0.9377609  0.44 0.455 0.4775 0.5266667 0.52875    39.43190 9.709012e-30 2.558368e-27 TRUE
# 36091      1 1031540 33 22 0.9550934  0.33 0.415 0.4475 0.4366667 0.41625    40.41913 6.398576e-32 2.402332e-18 TRUE
# 36167      1 1033999 39 22 0.9380956  0.41 0.425 0.4225 0.4566667 0.45125    39.79390 2.321236e-29 3.693890e-29 TRUE
# 37348      1 1065296 31 20 0.9111110  0.49 0.465 0.4600 0.4750000 0.51000    37.91483 7.670925e-29 9.673334e-19 TRUE
# 64452      1 1521805 48 24 0.9091562  0.53 0.540 0.5775 0.5400000 0.53625    37.38873 6.646896e-30 1.057444e-41 TRUE

names(result$fit)
# [1] "models"       "model_search" "fit.orig"     "fit"

result$fit$fit
#    chrom window_id start_index end_index start_pos   end_pos fit.k1 fit.k2     fit.f fit.isize     fit.llik
# 1      1       1:1           1     54732    943687 249136360      1      1 0.3619173     0.001 -322459.6812
# 2      2     2:548       54733    116870     50814 242852778      1      1 0.3619173     0.001 -366973.2160
# 3      3    3:1169      116871    173585    156233 197808975      1      1 0.3619173     0.001 -334137.2747
# 4      4    4:1736      173586    226485    367927 190789536      1      1 0.3619173     0.001 -311987.6236
# 5      5    5:2265      226486    276982    174940 180698588      1      1 0.3619173     0.001 -297948.0552
# 6      6    6:2770      276983    325496    231638 170800452      1      1 0.3619173     0.001 -286031.5297
# 7      7    7:3255      325497    365965     95280 158954519      1      1 0.3619173     0.001 -238868.5292
# 8      8    8:3660      365966    408725    313173 146173636      1      1 0.3619173     0.001 -252285.4733
# 9      9    9:4088      408726    438206    206838 140890997      1      1 0.3619173     0.001 -174087.5470
# 10    10   10:4383      438207    473436    159404 135433387      1      1 0.3619173     0.001 -207940.9852
# 11    11   11:4735      473437    508805    247630 134789668      1      1 0.3619173     0.001 -208765.5184
# 12    12   12:5089      508806    542562    269531 133501081      1      1 0.3619173     0.001 -198712.0743
# 13    13   13:5426      542563    549499  19509358  40757910      1      1 0.3619173     0.001  -40911.0420
# 14    13   13:5496      549500    549699  40758208  42306139      1      0 0.3619173     0.001   -1267.9631
# 15    13   13:5498      549700    549899  42307184  42708231      1      1 0.3619173     0.001   -1194.2028
# 16    13   13:5500      549900    549999  42711664  42886868      1      0 0.3619173     0.001    -628.9827
# 17    13   13:5501      550000    551199  42890109  46649721      1      1 0.3619173     0.001   -7030.2027
# 18    13   13:5513      551200    552299  46649723  50646218      1      0 0.3619173     0.001   -6739.4270
# 19    13   13:5524      552300    552499  50656582  51330918      0      0 0.3619173     0.001   -1116.2335
# 20    13   13:5526      552500    552699  51331121  51872361      1      0 0.3619173     0.001   -1249.8581
# 21    13   13:5528      552700    570664  51880895 115044332      1      1 0.3619173     0.001 -106131.0163
# 22    14   14:5707      570665    592761  20460865 107213845      1      1 0.3619173     0.001 -130438.8609
# 23    15   15:5928      592762    613627  20192951 102393537      1      1 0.3619173     0.001 -122951.4360
# 24    16   16:6137      613628    632327     94535  90115456      1      1 0.3619173     0.001 -109915.0982
# 25    17   17:6324      632328    647858     73263  81078768      1      1 0.3619173     0.001  -91272.3330
# 26    18   18:6479      647859    669309    132649  77967972      1      1 0.3619173     0.001 -126365.0313
# 27    19   19:6694      669310    678023    371967  59003830      1      1 0.3619173     0.001  -51133.5971
# 28    20   20:6781      678024    693002     61795  62884613      1      1 0.3619173     0.001  -88263.3816
# 29    21   21:6931      693003    702062  15452496  48052838      1      1 0.3619173     0.001  -53404.6820
# 30    22   22:7021      702063    707736  17284657  51151724      1      1 0.3619173     0.001  -33437.1102
```

## Future
* soma-cnv currently has a specific failure mode in which germ-line amplifications are mis-called as somatic CNV.  Future work will extend the error model to specifically include this possibility and reduce the incidence of germline-driven false positives.

