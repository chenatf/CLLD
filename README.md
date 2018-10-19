# CLLD (CRISPR large deletion Detector)
### version 1.0.0
## Description
This program is designed to detect germline de novo large deletions at the target sites using the NGS whole genome sequence data. 
## Authors
Rui Chen at the southern university of science and technology in the Jiankui He lab wrote the key source codes.
## Quick start
1. Installation.
```
git clone https://github.com/chenatf/CLLD.git
``` 
2.	Run the script.

* Mapping the raw sequence data to the reference
```
./cldd_mapping -r @RG\tID:id\tSM:sample\tLB:lib -f ref.fa -1 read1_lane1.fq,read1_lane2.fq -2 read2_lane1.fq,read2_lane2.fq -p sample -o output
```
* call germline de novo large deletion
```
./cldd_call -f ref.fa -s targets.bed -b parents1.bam,parents2.bam,offspring1.bam,offspring2.bam -p family -o output -a parents1,parents2 -g offspring1,offspring2
```

## Installation

This program is tested under Ubuntu 16.04 LTS.

1. Use the conda to easy install the packages.
```
conda install fastp
conda install bwa
conda install samtools
conda install mosdepth
conda install pandas
conda install pyvcf
conda install manta
```
2. download the speedseq from the github.
```
git clone --recursive https://github.com/hall-lab/speedseq
cd speedseq
make
```
3. download the program from the github.
```
git clone https://github.com/chenatf/CLLD.git
```

## Configuration

#### IO config
```
clld_io.config   record the input and output path of the program. Please do not change it if not necessary.
```
#### Software path config
```
clld_software.config   record the software path of the program. Please change it into your own path.
```

## Usage
CLLD is an analysis pipeline consists of two part:
1.	clld_mapping mapping the raw sequence data pair-end fastq files to the reference and get a duplicate-marked, sorted, indexed BAM file with speedseq align. It also calculate the mapping ratio and give the quality control and coverage analysis report.
2.	Clld_call read the family BAM files and call the large deletion in the genome with manta, then report the germline de novo large deletions happened in the target sites.



### clld_mapping
```
clld_mapping [options] -r <flag> -f <reference> -1 <pair1> -2 <pair2> -p <prefix> -o <output_dir> 
```
#### Required arguments
```
-r  STR  Read group header line such as "@RG\tID:id\tSM:sample\tLB:lib".
-f       Genome reference fasta file.
-1       Input paired-end fastq file(s) (comma separated).
         Example: -1 read1_lane1.fq,read1_lane2.fq,read1_lane3.fq.
-2       Input paired-end fastq file(s) (comma separated).
         Example: -2 read2_lane1.fq,read2_lane2.fq,read2_lane3.fq.
-p  STR  Output prefix
-o  STR  Output directory
 ```

#### Option arguments
```
-t  INT  Number of threads to use [default: 12]
-h       Show the usage messages.
```

#### Output
```
prefix.bam   The full, duplicate-marked, sorted BAM file for the library. This file may serve as input for clld_call.
prefix.html   Quality control report for fastq files by fastp.
prefix.dist.html   Coverage report for BAM files by mosdepth.
```

### clld_call
```
clld call [options] -f <reference> -s <bed> -b <bam> -p <prefix> -o <output_dir> -a <parents> -g <offspring>
```
#### Required arguments
```
-f      Genome reference fasta file.
-s      Target sites bed file. Three columns : chromosome start end.
-b      BAM file(s) (comma separated). example: -b in1.bam,in2.bam,in3.bam
-p      Output prefix
-o      Output directory
-a      Sample name of parents (comma separated), these should be consistent with the sample in the @RG.
	Example: parents1,parents2
-g      Sample name of offspring (comma separated), these should be consistent with the sample in the @RG.
	Example: offspring1,offspring2,offspring3
```

#### Option arguments
```
-t  INT  Number of threads to use [default: 24]
-h       Show the usage messages.
```
#### Output
```
Prefix.diploidSV.vcf   SVs and indels scored and genotyped under a diploid model for the set of samples in a joint diploid sample analysis by manta.
evidence.*.bam   Evidence reads of the candidate SVs identified from that input bam. Each read in an evidence bam keeps all information 
from the original bam by manta.
```
```
prefix.csv   A table report the germline de novo large deletions at the target sites.
```
Column | Description
------ | -----------
chromosome | chromosome of the large deletion
start	| Start site of the large deletion
end	| End site of the large deletion
length	| Length of the large deletion
Sample	| Sample name of the offspring
Target	| Target site of the large deletion
genotype	| Genotype of the large deletion

## Changelog

## Acknowledgements

## License and citation
This program is under the MIT license.
