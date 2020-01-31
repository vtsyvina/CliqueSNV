# CliqueSNV
## How to Run

Download jar from <a href="https://drive.google.com/open?id=1wAC57YARnVKi5uGAtPQwEHqAaPO_ru1m">here</a> (latest ver 1.4.11, January, 2020)

## Citation
Please cite preprint at BioRxiv: https://www.biorxiv.org/content/10.1101/264242v1

## Stable releases
1.4.1 - <a href="https://drive.google.com/file/d/18rVFT6fSSBdvHyPw0aibmlp3CfOsLV3K/view?usp=sharing">12 February 2018</a>

1.4.11 - <a href="https://drive.google.com/open?id=1Qwg6g9_lKlqY9zrH37f-F1QF32fWOmYt">2 January 2020</a>

## Parameters
There are several available parameters:
- ``-m`` mandatory parameter to specify method that you want to run. There are 4 posible values:
  - 'snv-pacbio' - run CliqueSNV with PacBio input
  - 'snv-illumina' - run CliqueSNV with Illumina input
  - 'snv-pacbio-vc' - run Variant Calling with CliqueSNV with PacBio input(.vcf file is output)
  - 'snv-illumina-vc' - run Variant Calling with CliqueSNV with Illumina input(.vcf file is output)
- ``-in`` input path. If not specified default sam files from ``data\PacBio_reads`` ro ``data\Illumina_reads`` folder will be used.
  It can be relative as well as absolute path. Illumina requires .sam file for both haplotyping and VC, PacBio can read
  .sam and .fas for haplotyping and .sam for VC. **Note!** if you use .fas for PacBio input, each read should be the same length.
  If it has offset it should be filled with 'N':

  ```
  >read1
  AACCTTGG
  >read2
  NACGTNNN
  ```
- ``-outDir`` output directory. `snv_output/` is default value
- ``-t`` - minimum threshold for O22 value. Default is 10 (only for Illumina reads)
- ``-tf`` - minimum threshold for 022 frequency relative to the reads' coverage. Default value os 0.000_333(10 out of 30 000)
- ``-log`` - some log data will be in console with this flag(no value needed)
- ``-cm`` - cliques merging algorithm. Two values: 'accurate', 'fast'. Default is 'accurate'. Accurate may lead to exponential explosion of cliques number.
That's why with large number of SNPs it may be useful to use fast polynomial algorithm with lower quality.
- ``-threads`` - number of threads for parallel tasks. By default program will use all available processors' cores.
- ``-os`` - output start position. If provided will cut the output from 0 to given position in haplotypes, variant calling
- ``-oe`` - output end position. If provided will cut the output from given position till the end in haplotypes, variant calling. 
For example, ``-os 100 -oe 700`` will output haplotypes only for positions [100, 700] or include SNPs in variant calling only inside this range
### t and tf parameters choice
These two parameters are significant, since they put a border in trade-off between precision and recall. By default they are set to detect very rare variants (<0.2%). If the data is known to be noisy(have a lot of read errors) and only variants with frequency >1% are of interest, then **-tf** should be around **0.01**, **-t** is optional and based on coverage.

### Usage example

``java -jar clique-snv.jar -m snv-pacbio``


``java -jar clique-snv.jar -m snv-illumina``(unzip sam file beforehand from 'data' folder)


``java -jar clique-snv.jar -m snv-illumina -in /path/to/data/r.sam -log``

``java -jar clique-snv.jar -m snv-illumina-vc -in /path/to/data/r.sam -outDir vcf_out/ -t 10 -tf 0.00034 -threads 8 -log``

### Example datasets
There are two example datasets:
- PacBio flu reads with 10 haplotypes
- Illumina flu reads with 2 haplotypes

``data/flu_ref.fasta`` contains those haplotypes as ground truth

How to run:

```java -jar clique-snv.jar -m snv-pacbio -log -in data\PacBio_reads\reads.sam```

```java -jar clique-snv.jar -m snv-illumina -in data\Illumina_reads\reads.sam```

## Output
**For default quasispecies problem** As output CliqueSNV has two files: human readable and fasta. Humanreadable file has the following form:
```
SNV got 10 haplotypes
[{
 snps=[],
 sourse clique='[],
 frequency= 0.5275826822460317,
 haplotype='GGAAAGAATAAAAGAACTAAGGAATCTAA...
'}, {
 snps=GT-TTATTAC[31, 265, 288, 396, 617, 747, 997, 1120, 1147, 2013],
 sourse clique='GTTATGTTAATTACC[31, 265, 267, 287, 288, 289, 396, 617, 747, 749, 997, 1120, 1147, 2013, 2014],
 frequency= 0.23675737005763478,
 haplotype='GGAAAGAATAAAAGAACTAAGGAATCTAATGG...
 '}, {
 ...
```
- 'snps' means alleles that differ from consensus in a certain found variant,
- 'source clique' is a service information specific to CliqueSNV method,
- 'frequency' - is an estimated frequency of a found variant in a range [0:1],
- 'haplotype' - is a full sequence for found variant

Fasta file will be:
```
>0_fr_0.5820184401895632
CCACAGCACGCAGATTGGTGGAATAAGGATGGTAAACATCCTTAGGCAGAACCC....

>1_fr_0.24979076133465727
CCACAGCACGCAGATTGGTGGAATAAGGATGGTAAACATCCTTAGGCAGAACCC...
 ...
```

Where name is just an index + haplotype frequency


**For Variant Calling** problem program produces standard VCF file. Standard is described <a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">here</a>

## Versions:
1.1.0 - add allele frequency for variant calling

1.2.0 - new cliques merging strategy; change true frequency estimator

1.3.0
- handle case when minor has frequency > 45% and can be actually major. Tool will automatically choose best option
- hence '-ref' argument is not maintained any more.
- performance improvement
- CliqueSNV will filter out haplotypes with frequencies < 8e-4 as false positives
- Fasta output for haplotyping

1.3.1 - parallel execution for Illumina input preprocessing

1.3.2
- optimized 'rotate minors with high frequency' step to make it in one iteration
- correctly handle single-paired reads

1.4.0
- new algorithm for merging cliques that allows to split one clique into several haplotypes considering all edges and 'no edges' between cliques
- improved performance for Bron-Kerbosch algorithm

1.4.1
- fix bug with clustering for reads with tie distance to several cliques

1.4.2
- add '-threads' parameter
- fix bugs for edge cases when no reads cover a position

1.4.3
- experimental implementation for imputation problem based on CliqueSNV added

1.4.4
- minor fix for edge case in clusters building

1.4.5
- further improvements of imputation algorithm (assumptions for read technology were removed)

1.4.6
- imputation method now can handle directory as input
- optimized multi-threading mechanism

1.4.7
- Fix bug when it is only one haplotype in the sample

1.4.8
- Handle consensus haplotype correctly; fix bug when the program finds only one haplotype

1.4.9
- Fix bug with exception for long reads in PacBio

1.4.10
- Add "-os", "-oe" parameters to control output range; small fixes

1.4.11
- Fix missing frequencies in output

## Any questions
With any questions. please, contact: vyacheslav.tsivina@gmail.com
