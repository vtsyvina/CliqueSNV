# CliqueSNV
## How to Run

Download jar from <a href="https://drive.google.com/file/d/1wAC57YARnVKi5uGAtPQwEHqAaPO_ru1m/view?usp=sharing">here</a> (latest ver 2.0.2, May 2021)

## How to build

`` mvn clean install``

It will create `clique-snv.jar` in the root folder

## Citation
Knyazev S, Tsyvina V, Shankar A, Melnyk A, Artyomenko A, Malygina T, Porozov YB, Campbell EM, Switzer WM, Skums P, Mangul S, Zelikovsky A. Accurate assembly of minority viral haplotypes from next-generation sequencing through efficient noise reduction. Nucleic Acids Res. 2021 Jul 2:gkab576. doi: 10.1093/nar/gkab576. Epub ahead of print. PMID: 34214168.

https://pubmed.ncbi.nlm.nih.gov/34214168/

## Stable releases
1.4.1 - <a href="https://drive.google.com/file/d/18rVFT6fSSBdvHyPw0aibmlp3CfOsLV3K/view?usp=sharing">12 February 2018</a>

1.4.11 - <a href="https://drive.google.com/open?id=1Qwg6g9_lKlqY9zrH37f-F1QF32fWOmYt">2 January 2020</a>

1.5.3 - <a href="https://drive.google.com/file/d/1MZ-J0O6yxSoF9OZxD6DT6ChazvOZsUa-"> 1 June 2020</a>

## Parameters
There are several available parameters:
- ``-m`` mandatory parameter to specify method that you want to run. There are 4 posible values:
  - 'snv-pacbio' - run CliqueSNV with PacBio input
  - 'snv-illumina' - run CliqueSNV with Illumina input
  - 'snv-pacbio-vc' - run Variant Calling with CliqueSNV with PacBio input(.vcf file is output)
  - 'snv-illumina-vc' - run Variant Calling with CliqueSNV with Illumina input(.vcf file is output)
  - 'consensus-illumina' - just a utility method to calculate consensus string base of input sam file for Illumina reads
  - 'consensus-pacbio' - just a utility method to calculate consensus string base of input sam file for Illumina reads
  
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

- ``-t`` - minimum threshold for O22 value. Default is 10
- ``-tf`` - minimum threshold for 022 frequency relative to the reads' coverage. Default value os 0.05. For more sensitive algorithm work decrease this parameter (may significantly increase runtime for diverse samples). **Note** Haplotypes with frequency <tf won't get into output
- ``-el`` - edges limit in SNPs graph. More edges will make the algorithm more sensitive. Too much links in the graph can lead to exponential explosion in time when calculating cliques. In version 1.5.0 only up to -el most frequent edges(SNPs supported by the most number of reads) are chosen. Default is 700(no limit in versions prior to 1.5.0) - should work good for the absolute majority of samples.  
- ``-sp`` - start position to search for SNPs (0-based, inclusive). If you hae a large genome, but is interested only in SNPS in specific region you can specify the range. It will make the program faster and more precise. 
But the output haplotypes still will be the whole genome(you can cut them using -os -oe parameters) - for example see below. Default value is 0
- ``-ep`` - end position to search for SNPs (0-based, inclusive). See -sp above
- ``-log`` - some log data will be in console with this flag(no value needed)
- ``-cm`` - cliques merging algorithm. Two values: 'accurate', 'fast'. Default is 'accurate'. Accurate may lead to exponential explosion of cliques number.
That's why with large number of SNPs it may be useful to use fast polynomial algorithm with lower quality.
- ``-ignoreDeletion`` - a flag to ignore deletions as potential SNVs for Illumina. For PacBio deletions are always ignored. 
- ``-threads`` - number of threads for parallel tasks. By default program will use all available processors' cores.
- ``-tl`` - the time limit in seconds. The program will stop after reaching this time limit and report just consensus haplotype. The json file will have "Timeout reached" under "error" in this case.
Default value is 10800 (3 hours). In most of the cases the program cannot finish because of too diverse sample with too much correlated SNPs (especially of long >10kb references).

### Output parameters
- ``-outDir`` output directory. `snv_output/` is default value
- ``-fdf`` fasta defline format. There are two options available(default is "short"):
    - ``short`` corresponds to ``>{id}_fr_{frequency}`` (e.g. ``>1_fr_0.5820184401895632``)
    - ``extended`` corresponds to ``>{sample_name}_{id}_{frequency}`` (e.g. ``>HIV_sample_1_0.92``). By default, the precision for 
    frequency is 2, it can be customized by adding the number to extended. For example, ``-fdf extended4`` will output ``>HIV_sample_1_0.9239``.
- ``-os`` - output start position(inclusive). If provided will cut the output from 0 to given position in haplotypes, variant calling
- ``-oe`` - output end position(exclusive). If provided will cut the output from given position till the end in haplotypes, variant calling. 
For example, ``-os 100 -oe 700`` will output haplotypes only for positions [100, 700) or include SNPs in variant calling only inside this range
- ``-rn`` - (no value needed) if present will create a file "<sample_name>_read_names.txt" with the reads assigned by the tool to all found haplotypes.
Format is the following:
```
><1_haplotype_name>
read1
read2
...
><2_haplotype_name>
read100
read102
...
```
**Note!** the same read can appear in different haplotypes if they are of the same distance to haplotypes. The tool assigns reads in several stages and the process is not trivial(for Illumina reads).
Thus, this file cannot be viewed as very precise as it is, but can be used for additional analysis.

### t and tf parameters choice
These two parameters are significant, since they put a border in trade-off between precision and recall. 
By default, they are set to detect moderate haplotypes (>5%). If it is know that data is not very noisy and variants with frequency >1% are of interest, then **-tf** should be around **0.01**, **-t** is optional and based on coverage.

### Difference between -sp -ep and -os -oe parameters
The (-sp,-ep) range doesn't affect the output, it affects only the region where the tool will work. And (-os, -oe) will cut the output without any affect on the program workflow.
Example:
if there are four haplotypes (and none of sp,ep,os,oe parameters specified):
```
AAAAAA
ACAAAC
AACCAA
AATGAA
```
(-sp=2, -ep=4). Output (second one disappear since SNPs are out of range and it won't be discovered):
```
AAAAAA
AACCAA
AATGAA
```

(-sp=2, -ep=4, -os=2, -oe=5):
```
AAA
ACC
TGA
```
(-os=2, -oe=5). It will only cut the haplotypes so you will see the duplicate here:
```
AAA
AAA
CCA
TGA
```
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

### Memory usage
From our experience the tool consumes around 10Gb(upper bound estimate) of RAM per 1,000,000 input reads(may vary based on a number of factors). To change standard JVM heap size limit specify -Xmx flag. Example with 50Gb:

```java -Xmx50G -jar clique-snv.jar -m snv-illumina -in data\Illumina_reads\reads.sam```

## Output
**For default quasispecies problem** As output CliqueSNV has two files: json and fasta. Json file has the info of used parameters, CliqueSNV version, found haplotypes:
```
{
  "version": "1.5.5",
  "settings": {
    "-m": "snv-pacbio",
    "-log": "true",
    "-in": "data\\PacBio_reads\\reads.sam",
    "-t": "10",
    "-rn": "true",
    "-tf": "0.0001"
  },
  "error": "none",
  "foundHaplotypes": 10,
  "haplotypes": [
    {
      "frequency": 0.5275842396392558,
      "name": "\u003e0_fr_0.5275842396392558",
      "snps": "[]",
      "sourceClique": "[]",
      "haplotype": "GGAAAGAATAAAAGAACTAAGGAATCTAA..."
    },
    {
      "frequency": 0.23674173677646565,
      "name": "\u003e1_fr_0.23674173677646565",
      "snps": "GT-TTATTAC[31, 265, 288, 396, 617, 747, 997, 1120, 1147, 2013]",
      "sourceClique": "GT-TTATTAC[31, 265, 288, 396, 617, 747, 997, 1120, 1147, 2013]",
      "haplotype": "GGAAAGAATAAAAGAACTAAGGAATCTAA..."
    },
 ...
```
- 'snps' means alleles that differ from consensus in a certain found variant,
- 'source clique' is a service information specific to CliqueSNV method,
- 'frequency' - is an estimated frequency of a found variant in a range [0:1],
- 'haplotype' - is a full sequence for found variant

Fasta file will be:
```
>1_fr_0.5820184401895632
CCACAGCACGCAGATTGGTGGAATAAGGATGGTAAACATCCTTAGGCAGAACCC....

>2_fr_0.24979076133465727
CCACAGCACGCAGATTGGTGGAATAAGGATGGTAAACATCCTTAGGCAGAACCC...
 ...
```

Where name is just an index + haplotype frequency


**For Variant Calling** problem program produces standard VCF file. Standard is described <a href="https://samtools.github.io/hts-specs/VCFv4.2.pdf">here</a>

## New in version 2.0.0 (sliding window approach)
CliqueSNV could perform poorly for samples with long references and significant diversity because of many linked SNPs for Illumina data.
To address it in version 2.0.0, we introduced a new approach to build haplotypes in windows, increasing the size of the window on each iteration.
The window size is equal to the sample's fragment length (usually 250 for single-end reads and around 500 for pair-end reads). 
The haplotypes reconstruction starts in the region of the highest coverage end extend the window to the left or right based on coverage.

This approach allows CliqueSNV to work on more diverse samples or samples with longer reference - 
we tested the new version on benchmarks with a reference length of 5000-10000, and CliqueSNV was able to reconstruct ground truth haplotypes successfully.


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

1.5.0
- Change of -tf parameter. Now default value is 5% and all haplotypes with frequency <tf value will be omitted
- Parallelization for several parts of algorithm that could take time for a big datasets(>2M reads) before. Should improve
the runtime on powerful machines
- Better handling of inputs with diverse quasispecies population(when pair-wise distances are >25-30 SNPs and many haplotypes within the sample)
- Output settings to console and to output file

1.5.1
- Start haplotypes from 1 in fasta output

1.5.2
- New -fdf parameter to change output fasta defline format

1.5.3
- Improved handling of SNPs graph for cliques search (should significantly improve performance on long references >5000 long)

1.5.4
- Better memory management, some performance improvements
- new "-rn" output parameter

1.5.5
- Change output from txt to json for ease of parsing
- "-tl"  parameter to set time limit

1.5.6
- "-sp", "-ep" parameters
- more error codes and messages in JSON output

1.5.7
- "-fc", "-ch" parameters to provide less conservative strategy as in versions 1.4.*
- put 'N' instead of '-' in haplotypes if the position is not covered by any read

2.0.0
- New "sliding window" approach in haplotypes reconstructing. More details above
- Small changes to handle deletion regions better and SNPs on the edges of reads
- "-fc", "-ch" - parameters not supported because of a new approach

## Any questions
With any questions. please, contact: v.tsyvina@gmail.com
