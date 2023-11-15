# iEVA

iEVA is a python command line tool for the Integration and Extraction of Variant Attributes.

> **Problem**

It is a good practice in DNAseq bioinformatics pipelines to include multiple variant callers to reach maximum sensitivity. However, merging information from all these variant callers in a harmonised and clear way is quite tricky.

Existing tools such as bcftools or GATK CombineVariants can merge multiple vcf files but only the annotation from one of them will be reported in the final merged file. Moreover, bioinformaticians and data scientists might be interested in getting additional information for each variant in the vcf files from different sources. For example, one might be interested in getting the GC content around the variant or the sequence conservation score. This information can only be obtained by parsing different datasets and the reference fasta file with the input vcf files and, eventually, the BAM file of the sample used to generate the vcf.

> **Solution**

iEVA is a python based command line tool that enables users to:

1. Merge multiple VCF files from different variant callers on both single- and multi-sample level.
2. Include external sequence-specific information including features from fasta reference file such as conservation score of the sequence, GC content, homopolymer content etc.
3. Generate a merged VCF file containing all the annotation from EACH input vcf file, conserving all the useful annotation derived from each variant caller, in addition to additional information explained above.

# Installation

**iEVA** is a command line tool developed in Python 2.7. To install iEVA, type in your terminal window:

    $ git clone https://github.com/Matteodigg/iEVA

Or download it from [here](https://github.com/Matteodigg/iEVA).

Enter the **iEVA** directory and execute:

    $ python setup.py install

Done! iEVA is installed in `/usr/local/bin`. To use iEVA, simply type:

    $ iEVA -I path/to/input.vcf -O path/to/output.vcf -Ref path/to/reference.fasta -[Arg1] ...

Detailed explanation of iEVA command line arguments can be found in the [iEVA arguments](#options) section.

> **See Also:**
>
> If a problem occurs installing python modules in the `REQUIREMENTS` file, please refer to the [module homepage](http://pypi.python.org/pypi) and try to install them separately using pypi:
>
>     $ pip install module_name

# Introduction

iEVA is a Python command line tool that expands the number of informative attributes that are not commonly provided by variant callers. It requires a `.vcf` file from variant calling and generates a new enriched vcf file as output.

iEVA attributes are annotated in the *INFO* and *FORMAT* fields of the vcf output file (see [vcf specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)). iEVA extracts information for the variant in the input vcf file from two different sources:

1. **Reference.fasta file**: Nucleotide sequence composition around a called variant is extracted from the fasta reference file. The annotations are reported in the vcf *INFO* field.
2. **BAM file**: Information about reads mapping variant position, for samples in the vcf input file, is extracted from the bam file. The annotations are reported in the vcf *FORMAT* field.

With iEVA, you can:
- Report the reference sequence quality in terms of G-C content, repeated sequences, pseudo nucleotide composition, and more.
- Report genotype information for samples in the vcf file. You can extract useful information for each sample, like fraction of unmapped reads, not paired reads, duplicate reads, etc.
- Report *allele-specific* information, evaluating, for example, the number of reads supporting reference or alternate allele given a specific mapping quality or base quality threshold. These attributes help identify any bias affecting called variants.

Main advantages of iEVA:
- The output file is in vcf format. iEVA simply adds the extracted attributes, maintaining original tags in the input vcf file.
- Extraction is independent of the variant caller used in variant calling analysis. For proper usage of iEVA, see [Usage](#how_to) and [Prerequisites](#Prereq).
- iEVA is not a variant caller, saving time. However, its computational time depends on the number of requested samples for extraction, target length, and mean target coverage.

See [iEVA arguments](#options) for details about iEVA options and the [Tutorial and examples](#ex) section.

One of the strengths of iEVA is its simplicity of integration into bioinformatics pipelines, as shown in this figure:
![workflow](/docs/source/Workflow-iEVA.PNG)


# Usage

iEVA requires three arguments to extract attributes from the reference.fasta file:

    $ iEVA -I path/to/input.vcf -Ref path/to/reference.fasta -O path/to/out.vcf -[opts]

An additional argument is required to extract attributes from the bam file:

    $ iEVA -I path/to/input.vcf -Ref path/to/reference.fasta -L path/to/Bam_list -O path/to/out.vcf -[opts]

The `Bam_list` file contains paths to each sample's bam file listed in the `--input` argument, one per line:

    path/to/Sample1.bam
    path/to/Sample2.bam
    ...

An index `.bai` file is required to be in the same path as the `.bam` file.

> **Note:**
>
> If you are running the `pyfasta` module on `reference.fasta` for the first time, it will take some time to write the index file.

# Prerequisites

For proper usage, iEVA meets the following prerequisites:

1. **Vcf normalization**:
    - Actual iEVA release does not extract any allele-specific attribute on multi-allelic sites with multiple values in vcf *ALT* field.
    - It is considered a good standard to *normalize* and *left-align* a vcf file after variant calling. Further details about *normalization* can be found [here](http://genome.sph.umich.edu/wiki/Variant_Normalization).

    You can normalize and split multi-allelic sites in a vcf file, before using iEVA, with [bcftools norm](https://samtools.github.io/bcftools/bcftools.html) using the following command line option:

        $ bcftools norm -m -both -f reference.fasta INPUT.vcf > OUTPUT.split-and-norm.vcf

2. **Bam header and sample name**:
    - To use *bam extraction*, the sample name in the vcf *genotype* field needs to be the same as the `RG:SM` tag in the bam file.
    - To modify the bam header sample name, you can use [Picard tool](https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) or [Samtools](http://www.htslib.org/doc/samtools.html).

> **Warning:**
>
> iEVA informs about extraction progress with the `--verbose` option. To use this argument, you have to **sort** the vcf file.

# Tutorial and Examples

This is the tutorial for iEVA. A detailed [argument description](#options) is available.

Starting from a `Test-norm.vcf` normalized file, a `ucsc.hg19.fasta` reference, and a `bamlist.txt`, let's take a look at how to use iEVA.

## Extracting Attributes from fasta Reference

Given an input vcf file from a variant caller containing these variants:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample-01	Sample-02
    chr1	55524197	.	A	G	.	.	AC=3;AF=0.75;AN=4;DP=509	GT:AD:DP:GQ	0/1:168,130:298:99	1/1:2,209:211:99
    chr2	179511767	.	GTA	G	.	.	AC=2;AF=0.5;AN=4;DP=104	GT:AD:DP:GQ	0/1:27,19:57:99	0/1:16,12:36:99
    chr3	12633424	.	A	AT	.	.	AC=2;AF=0.5;AN=4;DP=191	GT:AD:DP:GQ	0/1:49,20:95:99	0/1:39,22:83:99
    chr6	76550839	.	G	GT	.	.	AC=3;AF=0.75;AN=4;DP=251	GT:AD:DP:GQ	0/1:66,74:140:99	1/1:2,109:111:99
    chr9	137686987	.	A	T	.	.	AC=1;AF=0.25;AN=4;DP=249	GT:AD:DP:GQ	0/0:47,0:47:99	0/1:102,100:202:99
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226	GT:AD:DP:GQ	0/1:37,47:84:99	1/1:2,39:41:78
    chr18	28669664	.	CTTAA	C	.	.	AC=1;AF=0.25;AN=4;DP=84	GT:AD:DP:GQ	0/0:45,0:45:99	0/1:13,22:35:99

we verify if any variant falls in a repeated sequence, reporting its length and the single repeated unit for a window size of 300. In addition, we are going to check if variant falls in a lower case sequence for RepeatMasker:

    $ iEVA -v -I Test-norm.vcf -O iEVA-Test.vcf -Ref ucsc.hg19.fasta -SR -SRL -RM -SRU -WS 300

        Extracting attributes on: chr1

        Extracting attributes on: chr2

        Extracting attributes on: chr3

        Extracting attributes on: chr6

        Extracting attributes on: chr9

        Extracting attributes on: chr18

        Exatraction: Done

> note:
> 
> To extract only sequence features, argument `-L` for bam list is not required.

Output file `iEVA-Test.vcf` is now annotated as follow in *INFO* field:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=SR,Number=1,Type=String,Description="Variant falls in a repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">
    ##INFO=<ID=SRL,Number=1,Type=Integer,Description="Length of repeated sequence (expressed as number of nucleotides) for SR tag">
    ##INFO=<ID=SRU,Number=1,Type=String,Description="Simple repeated sequence unit composing repeated sequence (SR)">
    ##INFO=<ID=RM,Number=1,Type=Integer,Description="Variant falls in a repeated sequence according to RepeatMasker tool. True=1, False=0">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample-01	Sample-02
    chr1	55524197	.	A	G	.	.	AC=3;AF=0.75;AN=4;DP=509;SR=0;SRL=.;SRU=.;RM=0	GT:AD:DP:GQ	0/1:168,130:298:99	1/1:2,209:211:99
    chr2	179511767	.	GTA	G	.	.	AC=2;AF=0.5;AN=4;DP=104;SR=1;SRL=22;SRU=TA;RM=1	GT:AD:DP:GQ	0/1:27,19:57:99	0/1:16,12:36:99
    chr3	12633424	.	A	AT	.	.	AC=2;AF=0.5;AN=4;DP=191;SR=2;SRL=16;SRU=T;RM=0	GT:AD:DP:GQ	0/1:49,20:95:99	0/1:39,22:83:99
    chr6	76550839	.	G	GT	.	.	AC=3;AF=0.75;AN=4;DP=251;SR=0;SRL=.;SRU=.;RM=0	GT:AD:DP:GQ	0/1:66,74:140:99	1/1:2,109:111:99
    chr9	137686987	.	A	T	.	.	AC=1;AF=0.25;AN=4;DP=249;SR=0;SRL=.;SRU=.;RM=0	GT:AD:DP:GQ	0/0:47,0:47:99	0/1:102,100:202:99
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226;SR=0;SRL=.;SRU=.;RM=0	GT:AD:DP:GQ	0/1:37,47:84:99	1/1:2,39:41:78
    chr18	28669664	.	CTTAA	C	.	.	AC=1;AF=0.25;AN=4;DP=84;SR=0;SRL=.;SRU=.;RM=0	GT:AD:DP:GQ	0/0:45,0:45:99	0/1:13,22:35:99

You can observe that two variants have been found to be in a repeated sequence: a simple repeated sequence (`SR=1` composed by TA nucleotides repetitions for a total lenght `SRL=22`) and a homopolymer sequence (`SR=2` of nucleotide T with total length `SRL=16`), respectively.

Try to add other information about sequence nucleotide composition using *Pseudo Nucleotide Composition*, *GC* content and reporting *variant class*. Use as input the previous annotated output `iEVA-Test.vcf` by simply adding requested arguments:

    $ iEVA -v -I Test-norm.vcf -O iEVA-Test.vcf -Ref ucsc.hg19.fasta -SR -SRL -RM -PNC -GC -WS 300 -VC

Alternatively, you can use as input the previous annotated output `iEVA-Test.vcf` by simply adding requested arguments as follow:

    $ iEVA -v -I iEVA-Test.vcf -O iEVA-Test-PNC-GC-VC.vcf -Ref ucsc.hg19.fasta -PNC -GC -WS 300 -VC

        Extracting attributes on: chr1

        Extracting attributes on: chr2

        Extracting attributes on: chr3

        Extracting attributes on: chr6

        Extracting attributes on: chr9

        Extracting attributes on: chr18

        Exatraction: Done

In both cases, we have this output file:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=SR,Number=1,Type=String,Description="Variant falls in a repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">
    ##INFO=<ID=SRL,Number=1,Type=Integer,Description="Length of repeated sequence (expressed as number of nucleotides) for SR tag">
    ##INFO=<ID=PNC,Number=16,Type=Float,Description="Pseudo Nucleotide sequence Composition using Kmer size of 2. Reported as: AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT">
    ##INFO=<ID=RM,Number=1,Type=Integer,Description="Variant falls in a repeated sequence according to RepeatMasker tool. True=1, False=0">
    ##INFO=<ID=GC,Number=1,Type=Float,Description="Percentage of GC content in sequence">
    ##INFO=<ID=VC,Number=1,Type=String,Description="Annotated variant class: SNV=snv, Insertion=Ins, Deletion=Del, SequenceAlteration=Alt">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample-01	Sample-02
    chr1	55524197	.	A	G	.	.	AC=3;AF=0.75;AN=4;DP=509;SR=0;SRL=.;PNC=0.017,0.043,0.08,0.03,0.067,0.11,0.04,0.087,0.057,0.08,0.073,0.067,0.03,0.067,0.083,0.07;RM=0;GC=57.807;VC=snv	GT:AD:DP:GQ	0/1:168,130:298:99	1/1:2,209:211:99
    chr2	179511767	.	GTA	G	.	.	AC=2;AF=0.5;AN=4;DP=104;SR=1;SRL=22;PNC=0.137,0.04,0.067,0.12,0.05,0.013,0.0,0.073,0.047,0.017,0.03,0.04,0.133,0.067,0.037,0.13;RM=1;GC=26.91;VC=Del	GT:AD:DP:GQ	0/1:27,19:57:99	0/1:16,12:36:99
    chr3	12633424	.	A	AT	.	.	AC=2;AF=0.5;AN=4;DP=191;SR=2;SRL=16;PNC=0.077,0.05,0.073,0.073,0.1,0.053,0.003,0.073,0.043,0.047,0.043,0.047,0.053,0.083,0.06,0.12;RM=0;GC=41.196;VC=Ins	GT:AD:DP:GQ	0/1:49,20:95:99	0/1:39,22:83:99
    chr6	76550839	.	G	GT	.	.	AC=3;AF=0.75;AN=4;DP=251;SR=0;SRL=.;PNC=0.11,0.037,0.057,0.103,0.03,0.02,0.003,0.057,0.063,0.013,0.05,0.06,0.103,0.04,0.073,0.18;RM=0;GC=29.568;VC=Ins	GT:AD:DP:GQ	0/1:66,74:140:99	1/1:2,109:111:99
    chr9	137686987	.	A	T	.	.	AC=1;AF=0.25;AN=4;DP=249;SR=0;SRL=.;PNC=0.037,0.05,0.087,0.02,0.067,0.1,0.037,0.07,0.073,0.07,0.15,0.05,0.017,0.057,0.07,0.047;RM=0;GC=61.794;VC=snv	GT:AD:DP:GQ	0/0:47,0:47:99	0/1:102,100:202:99
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226;SR=0;SRL=.;PNC=0.053,0.047,0.067,0.017,0.063,0.16,0.023,0.093,0.063,0.06,0.11,0.047,0.007,0.073,0.077,0.04;RM=0;GC=61.794;VC=Del	GT:AD:DP:GQ	0/1:37,47:84:99	1/1:2,39:41:78
    chr18	28669664	.	CTTAA	C	.	.	AC=1;AF=0.25;AN=4;DP=84;SR=0;SRL=.;PNC=0.12,0.02,0.087,0.143,0.057,0.023,0.003,0.037,0.063,0.033,0.027,0.047,0.13,0.043,0.053,0.113;RM=0;GC=28.904;VC=Del	GT:AD:DP:GQ	0/0:45,0:45:99	0/1:13,22:35:99

## Extracting Sample level informations

As concern sample level extraction, we are going to test different attributes on genotype field ``Sample-01`` and ``Sample-02`` of vcf file.

Write a `Bamlist.txt` file with bam path of those samples, one for each raw:

    path/to/Sample-01.bam
    path/to/Sample-02.bam

iEVA allows to check the overall quality of a called variants using information stored in sample bam file. For example, try to extract some reads information for each sample like fraction of unmapped reads, not paired reads, not proper paired reads, mapping quality 0 reads, strand bias in reads orientation and duplicate reads with `-UnMap`, `-NP`, `-NPP`, `-MQ0`, `-SBR` and `-TDR` options:

    $ iEVA -v -I iEVA-Test.vcf -O iEVA-Test-02.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -UnMap -NP -NPP -TDR -MQ0 -SBR

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr6

    Extracting attributes on: chr9

    Extracting attributes on: chr18

    Exatraction: Done

Output file has, in addition to sequence information previously extracted, requested sample specific attributes in *FORMAT* vcf field:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=SR,Number=1,Type=String,Description="Variant falls in a repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">
    ##INFO=<ID=SRL,Number=1,Type=Integer,Description="Length of repeated sequence (expressed as number of nucleotides) for SR tag">
    ##INFO=<ID=RM,Number=1,Type=Integer,Description="Variant falls in a repeated sequence according to RepeatMasker tool. True=1, False=0">
    ##INFO=<ID=VC,Number=1,Type=String,Description="Annotated variant class: SNV=snv, Insertion=Ins, Deletion=Del, SequenceAlteration=Alt">
    ##FORMAT=<ID=SBR,Number=1,Type=Float,Description="Fisher exact test to detect strand bias (R1+,R1-,R2+,R2-)">
    ##FORMAT=<ID=UnMap,Number=1,Type=Float,Description="Fraction of unmapped reads">
    ##FORMAT=<ID=MQ0,Number=1,Type=Float,Description="Fraction of reads mapping position with Mapping Quaility=0">
    ##FORMAT=<ID=NP,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not paired">
    ##FORMAT=<ID=NPP,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not proper paired">
    ##FORMAT=<ID=TDR,Number=1,Type=Integer,Description="Fraction of total reads mapping position marked as duplicate">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample-01	Sample-02
    chr1	55524197	.	A	G	.	.	AC=3;AF=0.75;AN=4;DP=509;SR=0;SRL=.;RM=0;VC=snv	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR	0/1:168,130:298:99:0.7298:0:0:0:0:0.1497	1/1:2,209:211:99:1.0:0:0:0:0:0.1245
    chr2	179511767	.	GTA	G	.	.	AC=2;AF=0.5;AN=4;DP=104;SR=1;SRL=22;RM=1;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR	0/1:27,19:57:99:0.3988:0:0:0:0.0169:0.0781	0/1:16,12:36:99:0.5148:0:0:0:0:0.0513
    chr3	12633424	.	A	AT	.	.	AC=2;AF=0.5;AN=4;DP=191;SR=2;SRL=16;RM=0;VC=Ins	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR	0/1:49,20:95:99:0.5388:0:0:0:0:0.1471	0/1:39,22:83:99:0.5263:0:0:0:0:0.08
    chr6	76550839	.	G	GT	.	.	AC=3;AF=0.75;AN=4;DP=251;SR=0;SRL=.;RM=0;VC=Ins	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR	0/1:66,74:140:99:0.1552:0:0:0:0:0.1097	1/1:2,109:111:99:0.23:0:0:0:0.0472:0.1167
    chr9	137686987	.	A	T	.	.	AC=1;AF=0.25;AN=4;DP=249;SR=0;SRL=.;RM=0;VC=snv	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR	0/0:47,0:47:99:0.2271:0:0:0:0:0.1863	0/1:102,100:202:99:0.5754:0:0:0:0:0.2222
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226;SR=0;SRL=.;RM=0;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR	0/1:37,47:84:99:0.7298:0:0:0:0:0.0828	1/1:2,39:41:78:0.0956:0:0:0:0:0.0408
    chr18	28669664	.	CTTAA	C	.	.	AC=1;AF=0.25;AN=4;DP=84;SR=0;SRL=.;RM=0;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR	0/0:45,0:45:99:0.1448:0:0:0:0.0135:0.039	0/1:13,22:35:99:0.0456:0:0:0:0:0.0714

In the following example we are going to extract attributes on genotype field. For example, try to add iEVA read depth and iEVA allele depth with a base quality threshold of 20 (18 for InDels) and a mapping quality of 60 to extract the most informative reads. In addition, try to add mean *REF* q-score and mean *ALT* q-score:

    $ iEVA -v -I iEVA-Test-02.vcf -O iEVA-Test-03.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -iDP -iAD -iQR -iQA -SNVmbq 20 -INDELmbq 18 -SNVmmq 60 -INDELmmq 60

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr6

    Extracting attributes on: chr9

    Extracting attributes on: chr18

    Exatraction: Done

We get this output `iEVA-Test-03.vcf` file:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=SR,Number=1,Type=String,Description="Variant falls in a repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">
    ##INFO=<ID=SRL,Number=1,Type=Integer,Description="Length of repeated sequence (expressed as number of nucleotides) for SR tag">
    ##INFO=<ID=RM,Number=1,Type=Integer,Description="Variant falls in a repeated sequence according to RepeatMasker tool. True=1, False=0">
    ##INFO=<ID=VC,Number=1,Type=String,Description="Annotated variant class: SNV=snv, Insertion=Ins, Deletion=Del, SequenceAlteration=Alt">
    ##FORMAT=<ID=SBR,Number=1,Type=Float,Description="Fisher exact test to detect strand bias (R1+,R1-,R2+,R2-)">
    ##FORMAT=<ID=UnMap,Number=1,Type=Float,Description="Fraction of unmapped reads">
    ##FORMAT=<ID=MQ0,Number=1,Type=Float,Description="Fraction of reads mapping position with Mapping Quaility=0">
    ##FORMAT=<ID=NP,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not paired">
    ##FORMAT=<ID=NPP,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not proper paired">
    ##FORMAT=<ID=TDR,Number=1,Type=Integer,Description="Fraction of total reads mapping position marked as duplicate">
    ##FORMAT=<ID=iDP,Number=1,Type=Integer,Description="iEVA read depth. Only proper paired, proper mapped and not duplicate reads are included.">
    ##FORMAT=<ID=iAD,Number=R,Type=Integer,Description="Allelic depth reported by iEVA as Ref,Alt">
    ##FORMAT=<ID=iQR,Number=1,Type=Float,Description="Mean Q-score for REF allele">
    ##FORMAT=<ID=iQA,Number=1,Type=Float,Description="Mean Q-score for ALT allele">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample-01	Sample-02
    chr1	55524197	.	A	G	.	.	AC=3;AF=0.75;AN=4;DP=509;SR=0;SRL=.;RM=0;VC=snv	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA	0/1:168,130:298:99:0.7298:0:0:0:0:0.1497:300:168,132:29.76:31.77	1/1:2,209:211:99:1.0:0:0:0:0:0.1245:209:2,207:30.5:31.86
    chr2	179511767	.	GTA	G	.	.	AC=2;AF=0.5;AN=4;DP=104;SR=1;SRL=22;RM=1;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA	0/1:27,19:57:99:0.3988:0:0:0:0.0169:0.0781:49:30,19:30.45:30.89	0/1:16,12:36:99:0.5148:0:0:0:0:0.0513:28:16,12:30.69:30.63
    chr3	12633424	.	A	AT	.	.	AC=2;AF=0.5;AN=4;DP=191;SR=2;SRL=16;RM=0;VC=Ins	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA	0/1:49,20:95:99:0.5388:0:0:0:0:0.1471:82:60,22:30.0:29.14	0/1:39,22:83:99:0.5263:0:0:0:0:0.08:71:49,22:30.54:30.64
    chr6	76550839	.	G	GT	.	.	AC=3;AF=0.75;AN=4;DP=251;SR=0;SRL=.;RM=0;VC=Ins	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA	0/1:66,74:140:99:0.1552:0:0:0:0:0.1097:138:67,71:32.75:29.77	1/1:2,109:111:99:0.23:0:0:0:0.0472:0.1167:101:2,99:34.0:30.55
    chr9	137686987	.	A	T	.	.	AC=1;AF=0.25;AN=4;DP=249;SR=0;SRL=.;RM=0;VC=snv	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA	0/0:47,0:47:99:0.2271:0:0:0:0:0.1863:326:326,0:31.63:.	0/1:102,100:202:99:0.5754:0:0:0:0:0.2222:202:102,100:31.75:31.73
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226;SR=0;SRL=.;RM=0;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA	0/1:37,47:84:99:0.7298:0:0:0:0:0.0828:48:21,27:31.01:31.89	1/1:2,39:41:78:0.0956:0:0:0:0:0.0408:37:0,37:.:32.0
    chr18	28669664	.	CTTAA	C	.	.	AC=1;AF=0.25;AN=4;DP=84;SR=0;SRL=.;RM=0;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA	0/0:45,0:45:99:0.1448:0:0:0:0.0135:0.039:65:65,0:31.42:.	0/1:13,22:35:99:0.0456:0:0:0:0:0.0714:34:13,21:32.0:32.74

Now, we are going to look at some other interesting allele specific options.

For example, try to extract information about fraction of clipped reads for *REF* and *ALT* allele. Moreover, we want to check if duplicate reads supporting reference and alternate allele are well-balanced or bias affected. So, try to add `-iCR` and `-iCA` for clipped reads information and `-iDR`, `-iDA` and its difference `-iDDup` for duplicate reads with default values for mapping and base quality threshold.

> *Warning*
>
> Remember that optional arguments for base quality threshold `SNVmbq` and mapping quality threshold `SNVmmq` affect all the *allele-specific* arguments. In this case, default values will be used.

    $ iEVA -v -I iEVA-Test-03.vcf -O iEVA-Test-04.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -iCR -iCA -iDR -iDA -iDDup

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr6

    Extracting attributes on: chr9

    Extracting attributes on: chr18

    Exatraction: Done

This is the output `iEVA-Test-04.vcf` vcf file:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=SR,Number=1,Type=String,Description="Variant falls in a repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">
    ##INFO=<ID=SRL,Number=1,Type=Integer,Description="Length of repeated sequence (expressed as number of nucleotides) for SR tag">
    ##INFO=<ID=RM,Number=1,Type=Integer,Description="Variant falls in a repeated sequence according to RepeatMasker tool. True=1, False=0">
    ##INFO=<ID=VC,Number=1,Type=String,Description="Annotated variant class: SNV=snv, Insertion=Ins, Deletion=Del, SequenceAlteration=Alt">
    ##FORMAT=<ID=SBR,Number=1,Type=Float,Description="Fisher exact test to detect strand bias (R1+,R1-,R2+,R2-)">
    ##FORMAT=<ID=UnMap,Number=1,Type=Float,Description="Fraction of unmapped reads">
    ##FORMAT=<ID=MQ0,Number=1,Type=Float,Description="Fraction of reads mapping position with Mapping Quaility=0">
    ##FORMAT=<ID=NP,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not paired">
    ##FORMAT=<ID=NPP,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as not proper paired">
    ##FORMAT=<ID=TDR,Number=1,Type=Integer,Description="Fraction of total reads mapping position marked as duplicate">
    ##FORMAT=<ID=iDP,Number=1,Type=Integer,Description="iEVA read depth. Only proper paired, proper mapped and not duplicate reads are included.">
    ##FORMAT=<ID=iAD,Number=R,Type=Integer,Description="Allelic depth reported by iEVA as Ref,Alt">
    ##FORMAT=<ID=iQR,Number=1,Type=Float,Description="Mean Q-score for REF allele">
    ##FORMAT=<ID=iQA,Number=1,Type=Float,Description="Mean Q-score for ALT allele">
    ##FORMAT=<ID=iDR,Number=1,Type=Float,Description="Fraction of duplicate reads mapping REF allele">
    ##FORMAT=<ID=iDA,Number=1,Type=Float,Description="Fraction of duplicate reads mapping ALT allele">
    ##FORMAT=<ID=iDDup,Number=1,Type=Float,Description="Difference between fraction of duplicate reads for REF and ALT alleles (DupREF-DupALT)">
    ##FORMAT=<ID=iClipRef,Number=1,Type=Float,Description="Fraction of clipped reads supporting REF">
    ##FORMAT=<ID=iClipAlt,Number=1,Type=Float,Description="Fraction of clipped reads supporting ALT">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample-01	Sample-02
    chr1	55524197	.	A	G	.	.	AC=3;AF=0.75;AN=4;DP=509;SR=0;SRL=.;RM=0;VC=snv	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA:iDR:iDA:iDDup:iCR:iCA	0/1:168,130:298:99:0.7298:0:0:0:0:0.1497:300:168,132:29.76:31.77:0.1508:0.1484:0.0024:0.0059:0.0152	1/1:2,209:211:99:1.0:0:0:0:0:0.1245:209:2,207:30.5:31.86:0:0.1255:-0.1255:0:0.0144
    chr2	179511767	.	GTA	G	.	.	AC=2;AF=0.5;AN=4;DP=104;SR=1;SRL=22;RM=1;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA:iDR:iDA:iDDup:iCR:iCA	0/1:27,19:57:99:0.3988:0:0:0:0.0169:0.0781:49:30,19:30.45:30.89:0.0606:0.0952:-0.0346:0:0	0/1:16,12:36:99:0.5148:0:0:0:0:0.0513:28:16,12:30.69:30.63:0.1111:0:0.1111:0:0
    chr3	12633424	.	A	AT	.	.	AC=2;AF=0.5;AN=4;DP=191;SR=2;SRL=16;RM=0;VC=Ins	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA:iDR:iDA:iDDup:iCR:iCA	0/1:49,20:95:99:0.5388:0:0:0:0:0.1471:82:60,22:30.0:29.14:0.1389:0.08:0.0589:0.0968:0.087	0/1:39,22:83:99:0.5263:0:0:0:0:0.08:71:49,22:30.54:30.64:0.0755:0.12:-0.0445:0.102:0
    chr6	76550839	.	G	GT	.	.	AC=3;AF=0.75;AN=4;DP=251;SR=0;SRL=.;RM=0;VC=Ins	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA:iDR:iDA:iDDup:iCR:iCA	0/1:66,74:140:99:0.1552:0:0:0:0:0.1097:138:67,71:32.75:29.77:0.141:0.0779:0.0631:0:0	1/1:2,109:111:99:0.23:0:0:0:0.0472:0.1167:101:2,99:34.0:30.55:0:0.1161:-0.1161:0:0
    chr9	137686987	.	A	T	.	.	AC=1;AF=0.25;AN=4;DP=249;SR=0;SRL=.;RM=0;VC=snv	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA:iDR:iDA:iDDup:iCR:iCA	0/0:47,0:47:99:0.2271:0:0:0:0:0.1863:326:326,0:31.63:.:0.1867:.:.:0.003:.	0/1:102,100:202:99:0.5754:0:0:0:0:0.2222:202:102,100:31.75:31.73:0.2031:0.2424:-0.0393:0:0
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226;SR=0;SRL=.;RM=0;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA:iDR:iDA:iDDup:iCR:iCA	0/1:37,47:84:99:0.7298:0:0:0:0:0.0828:48:21,27:31.01:31.89:0.0435:0.0882:-0.0447:0:0	1/1:2,39:41:78:0.0956:0:0:0:0:0.0408:37:0,37:.:32.0:.:0.025:.:.:0.1795
    chr18	28669664	.	CTTAA	C	.	.	AC=1;AF=0.25;AN=4;DP=84;SR=0;SRL=.;RM=0;VC=Del	GT:AD:DP:GQ:SBR:UnMap:MQ0:NP:NPP:TDR:iDP:iAD:iQR:iQA:iDR:iDA:iDDup:iCR:iCA	0/0:45,0:45:99:0.1448:0:0:0:0.0135:0.039:65:65,0:31.42:.:0.0441:.:.:0:.	0/1:13,22:35:99:0.0456:0:0:0:0:0.0714:34:13,21:32.0:32.74:0.1333:0.0455:0.0878:0:0

Finally, try to extract attributes both from reference fasta file and bam file only for `Sample-01`. In this case, simply remove `Sample-02` bam file from `Bamlist.txt`. Missing values for `Sample-02` are indicated by a `.` for all requested optional arguments.

    $ iEVA -v -I Test-norm.vcf -O iEVA-Test-Only-Sample-01.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -SR -SRL -RM -GC -WS 300 -AS -UnMap -SA -DDup -iDP -iAD -iQR -iQA -iRMQ -iAMQ

    Sample 'Sample-02' will not be annotated. Missing bam file.

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr6

    Extracting attributes on: chr9

    Extracting attributes on: chr18

    Exatraction: Done

We get the following result:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##INFO=<ID=SR,Number=1,Type=String,Description="Variant falls in a repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">
    ##INFO=<ID=SRL,Number=1,Type=Integer,Description="Length of repeated sequence (expressed as number of nucleotides) for SR tag">
    ##INFO=<ID=RM,Number=1,Type=Integer,Description="Variant falls in a repeated sequence according to RepeatMasker tool. True=1, False=0">
    ##INFO=<ID=GC,Number=1,Type=Float,Description="Percentage of GC content in sequence">
    ##FORMAT=<ID=UnMap,Number=1,Type=Float,Description="Fraction of unmapped reads">
    ##FORMAT=<ID=SA,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as supplementary alignment">
    ##FORMAT=<ID=AS,Number=1,Type=Float,Description="Reads mean alignment score">
    ##FORMAT=<ID=iDDup,Number=1,Type=Float,Description="Difference between fraction of duplicate reads for REF and ALT alleles (DupREF-DupALT)">
    ##FORMAT=<ID=iDP,Number=1,Type=Integer,Description="iEVA read depth. Only proper paired, proper mapped and not duplicate reads are included.">
    ##FORMAT=<ID=iAD,Number=R,Type=Integer,Description="Allelic depth reported by iEVA as Ref,Alt">
    ##FORMAT=<ID=iQR,Number=1,Type=Float,Description="Mean Q-score for REF allele">
    ##FORMAT=<ID=iQA,Number=1,Type=Float,Description="Mean Q-score for ALT allele">
    ##FORMAT=<ID=iRMQ,Number=1,Type=Float,Description="Mean mapping quality score for reads supporting REF allele">
    ##FORMAT=<ID=iAMQ,Number=1,Type=Float,Description="Mean mapping quality score for reads supporting ALT allele">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample-01	Sample-02
    chr1	55524197	.	A	G	8369.19	.	AC=3;AF=0.75;AN=4;DP=509;SR=0;SRL=.;RM=0;GC=57.807	GT:AD:DP:GQ:UnMap:SA:AS:iDDup:iDP:iAD:iQR:iQA:iRMQ:iAMQ	0/1:168,130:298:99:0:0:127.502:0.0024:301:169,132:29.67:31.77:60.0:60.0	1/1:2,209:211:99:.:.:.:.:.:.:.:.:.:.
    chr2	179511767	.	GTA	G	763.4	.	AC=2;AF=0.5;AN=4;DP=104;SR=1;SRL=22;RM=1;GC=26.91	GT:AD:DP:GQ:UnMap:SA:AS:iDDup:iDP:iAD:iQR:iQA:iRMQ:iAMQ	0/1:27,19:57:99:0:0:141.441:-0.0346:50:31,19:30.27:30.89:60.0:62.6316	0/1:16,12:36:99:.:.:.:.:.:.:.:.:.:.
    chr3	12633424	.	A	AT	553.4	.	AC=2;AF=0.5;AN=4;DP=191;SR=2;SRL=16;RM=0;GC=41.196	GT:AD:DP:GQ:UnMap:SA:AS:iDDup:iDP:iAD:iQR:iQA:iRMQ:iAMQ	0/1:49,20:95:99:0:0:133.871:0.0589:85:62,23:29.44:28.3:60.0:60.0	0/1:39,22:83:99:.:.:.:.:.:.:.:.:.:.
    chr6	76550839	.	G	GT	6860.15	.	AC=3;AF=0.75;AN=4;DP=251;SR=0;SRL=.;RM=0;GC=29.568	GT:AD:DP:GQ:UnMap:SA:AS:iDDup:iDP:iAD:iQR:iQA:iRMQ:iAMQ	0/1:66,74:140:99:0:0:142.928:0.0631:138:67,71:32.75:29.77:60.0:60.0	1/1:2,109:111:99:.:.:.:.:.:.:.:.:.:.
    chr9	137686987	.	A	T	2070.19	.	AC=1;AF=0.25;AN=4;DP=249;SR=0;SRL=.;RM=0;GC=61.794	GT:AD:DP:GQ:UnMap:SA:AS:iDDup:iDP:iAD:iQR:iQA:iRMQ:iAMQ	0/0:47,0:47:99:0:0:137.367:.:331:331,0:31.39:.:60.0:.	0/1:102,100:202:99:.:.:.:.:.:.:.:.:.:.
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226;SR=0;SRL=.;RM=0;GC=61.794	GT:AD:DP:GQ:UnMap:SA:AS:iDDup:iDP:iAD:iQR:iQA:iRMQ:iAMQ	0/1:37,47:84:99:0:0:126.602:-0.0447:53:22,31:30.89:30.5:60.0:60.0	1/1:2,39:41:78:.:.:.:.:.:.:.:.:.:.
    chr18	28669664	.	CTTAA	C	814.15	.	AC=1;AF=0.25;AN=4;DP=84;SR=0;SRL=.;RM=0;GC=28.904	GT:AD:DP:GQ:UnMap:SA:AS:iDDup:iDP:iAD:iQR:iQA:iRMQ:iAMQ	0/0:45,0:45:99:0:0:147.581:.:65:65,0:31.42:.:60.0:.	0/1:13,22:35:99:.:.:.:.:.:.:.:.:.:.

A really important feature of iEVA is its flexibility and usage on more variant callers. For example, try to merge results coming out from different variant caller, obtaining a vcf file like this with different variant caller-specific annotations:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">
    ##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
    ##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
    ##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
    ##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
    ##INFO=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">
    ##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 15">
    ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
    ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
    ##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
    ##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
    ##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
    ##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
    ##INFO=<ID=ADP,Number=1,Type=Integer,Description="Average per-sample depth of bases with Phred score >= 15">
    ##INFO=<ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">
    ##INFO=<ID=HET,Number=1,Type=Integer,Description="Number of samples called heterozygous-variant">
    ##INFO=<ID=HOM,Number=1,Type=Integer,Description="Number of samples called homozygous-variant">
    ##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of samples not called">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample-01	Sample-02
    chr1	55524197	.	A	G	.	.	AC=3;AF=0.75;AN=4;DP=509	GT:AD:DP:GQ	0/1:168,130:298:99	1/1:2,209:211:99
    chr1	55524197	.	A	G	10204	.	AC=3;AF=0.75;AN=4;AO=341;DP=512;QA=12500;QR=5633;TYPE=snp	GT:GQ:DP:AD:RO:QR:AO:QA	0/1:135.817:301:169,132:169:5566:132:4859	1/1:135.817:211:2,209:2:67:209:7641
    chr1	55524197	.	A	G	.	.	ADP=148;WT=0;HET=1;HOM=1	GT:GQ:DP:RD:AD:FREQ:PVAL:RBQ:ABQ	0/1:255:181:100:81:44,75%:3,5278E-30:55:60	1/1:255:116:1:115:99,14%:3,2397E-67:67:66
    chr2	179511767	.	GTA	G	.	.	AC=2;AF=0.5;AN=4;DP=104	GT:AD:DP:GQ	0/1:27,19:57:99	0/1:16,12:36:99
    chr2	179511767	.	GTA	G	0	.	AC=0;AF=0;AN=4;AO=32;DP=93;QA=995;QR=1410;TYPE=del	GT:GQ:DP:AD:RO:QR:AO:QA	0/0:160.002:56:27,19:27:885:19:606	0/0:160.002:37:16,13:16:525:13:389
    chr2	179511767	.	GTA	G	.	.	ADP=33;WT=0;HET=2;HOM=0	GT:GQ:DP:RD:AD:FREQ:PVAL:RBQ:ABQ	0/1:45:43:22:13:30,23%:3,1103E-5:47:53	0/1:24:24:11:7:29,17%:3,8123E-3:53:64
    chr3	12633424	.	A	AT	.	.	AC=2;AF=0.5;AN=4;DP=191	GT:AD:DP:GQ	0/1:49,20:95:99	0/1:39,22:83:99
    chr3	12633424	.	A	AT	704	.	AC=2;AF=0.5;AN=4;DP=194;QA=1402;QR=3200;TYPE=ins	GT:GQ:DP:AD:RO:QR:AO:QA	0/1:136.293:108:56,24:56:1821:24:699:	0/1:136.293:86:42,22:42:1379:22:703:
    chr3	12633424	.	A	AT	.	.	ADP=71;WT=0;HET=2;HOM=0	GT:GQ:DP:RD:AD:FREQ:PVAL:RBQ:ABQ	0/1:49:86:51:15:17,24%:1,2459E-5:44:45	0/1:47:56:29:14:25%:1,727E-5:57:51
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226	GT:AD:DP:GQ	0/1:37,47:84:99	1/1:2,39:41:78
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	795	.	AC=2;AF=0.5;AN=4;DP=120;QA=1612;QR=787;TYPE=del	GT:GQ:DP:AD:RO:QR:AO:QA	0/1:160.002:81:23,31:23:787:31:666	0/1:160.002:39:0,39:0:0:39:946:
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	ADP=69;WT=0;HET=2;HOM=0	GT:GQ:DP:RD:AD:FREQ:PVAL:RBQ:ABQ	0/1:41:80:67:13:16,25%:7,1891E-5:52:62	0/1:76:58:36:22:37,93%:2,025E-8:56:57

To retrieve sequence attributes and genotype information for both `Sample-01` and `Sample-02` use iEVA as showed in previous examples:

    $ iEVA -v -I merged-norm.vcf -O iEVA-merged.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -SR -SRL -RM -GC -WS 300 -AS -Unmap -SA -iDDup -MQ0 -iDP -iAD -iQR -iQA

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr9

    Exatraction: Done

As you can see, resulting file reports the input vcf file with iEVA attributes extracted for all variants, independently from the variant caller used:

    ##fileformat=VCFv4.2
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
    ##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">
    ##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Sum of quality of the alternate observations">
    ##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Sum of quality of the reference observations">
    ##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
    ##INFO=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">
    ##INFO=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">
    ##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Quality Read Depth of bases with Phred score >= 15">
    ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
    ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
    ##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
    ##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
    ##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description="Average quality of reference-supporting bases (qual1)">
    ##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description="Average quality of variant-supporting bases (qual2)">
    ##INFO=<ID=ADP,Number=1,Type=Integer,Description="Average per-sample depth of bases with Phred score >= 15">
    ##INFO=<ID=WT,Number=1,Type=Integer,Description="Number of samples called reference (wild-type)">
    ##INFO=<ID=HET,Number=1,Type=Integer,Description="Number of samples called heterozygous-variant">
    ##INFO=<ID=HOM,Number=1,Type=Integer,Description="Number of samples called homozygous-variant">
    ##INFO=<ID=NC,Number=1,Type=Integer,Description="Number of samples not called">
    ##INFO=<ID=SR,Number=1,Type=String,Description="Variant falls in a repeated sequence. None=0, SimpleRepeat=1, Homopolymer=2.">
    ##INFO=<ID=SRL,Number=1,Type=Integer,Description="Length of repeated sequence (expressed as number of nucleotides) for SR tag">
    ##INFO=<ID=RM,Number=1,Type=Integer,Description="Variant falls in a repeated sequence according to RepeatMasker tool. True=1, False=0">
    ##INFO=<ID=GC,Number=1,Type=Float,Description="Percentage of GC content in sequence">
    ##FORMAT=<ID=UnMap,Number=1,Type=Float,Description="Fraction of unmapped reads">
    ##FORMAT=<ID=MQ0,Number=1,Type=Float,Description="Fraction of reads mapping position with Mapping Quaility=0">
    ##FORMAT=<ID=SA,Number=1,Type=Float,Description="Fraction of reads mapping position flagged as supplementary alignment">
    ##FORMAT=<ID=AS,Number=1,Type=Float,Description="Reads mean alignment score">
    ##FORMAT=<ID=iDDup,Number=1,Type=Float,Description="Difference between fraction of duplicate reads for REF and ALT alleles (DupREF-DupALT)">
    ##FORMAT=<ID=iDP,Number=1,Type=Integer,Description="iEVA read depth. Only proper paired, proper mapped and not duplicate reads are included.">
    ##FORMAT=<ID=iAD,Number=R,Type=Integer,Description="Allelic depth reported by iEVA as Ref,Alt">
    ##FORMAT=<ID=iQR,Number=1,Type=Float,Description="Mean Q-score for REF allele">
    ##FORMAT=<ID=iQA,Number=1,Type=Float,Description="Mean Q-score for ALT allele">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	20161125_01_Cardio_Sort	20161125_02_Cardio_Sort
    chr1	55524197	.	A	G	.	.	AC=3;AF=0.75;AN=4;DP=509;SR=0;SRL=.;RM=0;GC=57.807	GT:AD:DP:GQ:UnMap:MQ0:SA:AS:iDDup:iDP:iAD:iQR:iQA	0/1:168,130:298:99:0:0:0:127.502:0.0024:301:169,132:29.67:31.77	1/1:2,209:211:99:0:0:0:107.332:-0.1255:211:2,209:30.5:31.7
    chr2	179511767	.	GTA	G	0	.	AC=0;AF=0;AN=4;AO=32;DP=93;QA=995;QR=1410;TYPE=del;SR=1;SRL=22;RM=1;GC=26.91	GT:GQ:DP:AD:RO:QR:AO:QA:UnMap:MQ0:SA:AS:iDDup:iDP:iAD:iQR:iQA	0/0:160.002:56:27,19:27:885:19:606:0:0:0:141.441:-0.0346:50:31,19:30.27:30.89	0/0:160.002:37:16,13:16:525:13:389:0:0:0:137.378:0.1111:29:16,13:30.69:30.0
    chr3	12633424	.	A	AT	.	.	ADP=71;WT=0;HET=2;HOM=0;SR=2;SRL=16;RM=0;GC=41.196	GT:GQ:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:UnMap:MQ0:SA:AS:iDDup:iDP:iAD:iQR:iQA	0/1:49:86:51:15:17,24%:1,2459E-5:44:45:0:0:0:133.871:0.0589:85:62,23:29.44:28.3	0/1:47:56:29:14:25%:1,727E-5:57:51:0:0:0:125.902:-0.0445:71:49,22:30.54:30.64
    chr9	137688780	.	GCCTGTCCCCTCCAAAACCCA	G	.	.	AC=3;AF=0.75;AN=4;DP=226;SR=0;SRL=.;RM=0;GC=61.794	GT:AD:DP:GQ:UnMap:MQ0:SA:AS:iDDup:iDP:iAD:iQR:iQA	0/1:37,47:84:99:0:0:0:126.602:-0.0447:53:22,31:30.89:30.5	1/1:2,39:41:78:0:0:0:123.053:.:39:0,39:.:31.46

# iEVA arguments

A full list of iEVA tool arguments with a short description can be accessed by typing in your command line `-h` [- -help] option:

    $ iEVA -h [--help]

    usage: iEVA [-h] [-v] -I path/to/input -Ref path/to/fasta -O path/to/outfile
                [-L path/to/bamlist] [-WS [20-600]] [-SR]
                [-SRList path/to/SampleRepeatList.txt] [-SRL] [-SRU] [-PNC] [-RM]
                [-GC] [-VC] [-SBR] [-UnMap] [-MMQ] [-MQ0] [-NPA] [-SA] [-NP]
                [-NPP] [-AS] [-TDR] [-iNDR] [-iNDA] [-iDR] [-iDA] [-iDDup] [-TDP]
                [-SNVmbq [0-66]] [-SNVmmq [0-60]] [-INDELmbq [0-66]]
                [-INDELmmq [0-60]] [-iDP] [-iAD] [-iRR] [-iRA] [-iQR] [-iQA]
                [-iRMQ] [-iAMQ] [-iNCR] [-iNCA] [-iCR] [-iCA]
    
    iEVA: identification and Extraction of Variant Attributes.
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         Increase output verbosity. Enable iEVA extraction
                            progress information
    
      Required arguments
    
      -I path/to/input, --input path/to/input
                            Input file in vcf format.
      -Ref path/to/fasta, --reference path/to/fasta
                            Reference file in fasta format.
      -O path/to/outfile, --outfile path/to/outfile
                            Output vcf file.
      -L path/to/bamlist, --list path/to/bamlist
                            REQUIRED ONLY TO EXTRACT FEATURES FROM BAM FILE. Bam
                            list file, including path, one for each raw. Sample
                            name in vcf genotype field has to be the same of Read
                            Group Sample Name (SM) tag in bam file header. Bam
                            Index file required in same path of bam file.
    
      Reference.fasta extraction arguments
    
      -WS [20-600], --WindowSize [20-600]
                            Set window size for -SR -SRL -GC and -PNC command.
                            Default=40, Min=20, Max=600
      -SR, --SimpleRepeat   Enable Simple Repeat Finder. Check if a variant fall
                            in a Simple Repeated Sequence. 0 for None, 1 for
                            Simple Repeat Sequence, 2 for Homopolymer.
      -SRList path/to/SampleRepeatList.txt, --SimpleRepeatList path/to/SampleRepeatList.txt
                            path to list containing simple repeated sequence (one
                            for each raw) to extract with iEVA. If enabled, only
                            reported sequences in file will be extracted. -SR
                            option is required.
      -SRL, --SimpleRepeatLength
                            Length of repeated sequence (expressed as number of
                            nucleotides) in -SR option
      -SRU, --SimpleRepeatUnit
                            Repeated sequence unit composing simple repeated
                            sequence
      -PNC, --PseudoNucleotidesComposition
                            Describe nucleotide sequence using Pseudo Nucleotide
                            Composition with Kmer size of 2. Values reported as:
                            AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT
      -RM, --RepeatMasker   Check from fasta reference if a variant falls in a
                            sequence flagged as repeated sequence by RepeatMasker
                            tool.
      -GC, --gcContent      Percentage of sequence GC content.
      -VC, --VariantClass   Annotated variant class: SNV (snv), Insertion (Ins),
                            Deletion (Del), Sequence Alteration (Alt)
    
      Bam extraction arguments
    
      -SBR, --StrandBiasReads
                            Fisher's exact test based on read orientation
                            (R1+,R1-,R2+,R2-) to detect strand bias. Only for
                            paired-end sequencing experiments.
      -UnMap, --UnMappedReads
                            Fraction of unmapped reads for variant position.
      -MMQ, --MeanMappingQuality
                            Mean mapping quality for reads mapping variant
                            position
      -MQ0, --MappingQualityZero
                            Fraction of reads with Mapping Quaility=0.
      -NPA, --NotPrimaryAlignment
                            Fraction of reads mapping position flagged as not
                            primary alignment.
      -SA, --SupplementaryAlignment
                            Fraction of reads mapping position flagged as
                            supplementary alignment.
      -NP, --NotPairedReads
                            Fraction of reads mapping position flagged as not
                            paired.
      -NPP, --NotProperPairedReads
                            Fraction of reads mapping position flagged as not
                            proper paired.
      -AS, --AlignmentScore
                            Mean Alignment Score of reads mapping position.
      -TDR, --TotalDupReads
                            Fraction of total reads mapping position marked as
                            duplicate.
    
      Genotype extraction arguments
    
      -iNDR, --NumberReadDupRef
                            Number of reads mapping position marked as duplicate
                            on REF allele.
      -iNDA, --NumberReadDupAlt
                            Number of reads mapping position marked as duplicate
                            on ALT allele.
      -iDR, --DuplicateReference
                            Fraction of reads mapping position marked as duplicate
                            on REF allele.
      -iDA, --DuplicateAlternate
                            Fraction of reads mapping position marked as duplicate
                            on ALT allele.
      -iDDup, --DeltaDuplicate
                            Difference between fraction duplicate reads in REF and
                            ALT alleles (REF-ALT). Negative values show preference
                            in duplicate for REF allele, positive for ALT.
      -TDP, --TotalDPUnfilter
                            Total read depth. No filter applied, include duplicate
                            reads. Look at -iDP option for filtered iDP.
      -SNVmbq [0-66], --SNVMinBaseQuality [0-66]
                            Minimum Base Quality threshold for base supporting SNV
                            position. Used on allele specific annotation.
                            Default=12
      -SNVmmq [0-60], --SNVMinMappingQuality [0-60]
                            Minimum Mapping Quality threshold for reads supporting
                            SNV position. Used on allele specific annotation.
                            Default=30
      -INDELmbq [0-66], --IndelMinBaseQuality [0-66]
                            Minimum Base Quality threshold for base supporting
                            InDel position. Used on allele specific annotation.
                            Default=10
      -INDELmmq [0-60], --IndelMinMappingQuality [0-60]
                            Minimum Mapping Quality threshold for reads supporting
                            InDel position. Used on allele specific annotation.
                            Default=20
      -iDP, --iEvaDepth     iEVA read depth for variant position. Only proper
                            paired, proper mapped and not duplicate reads are
                            included.
      -iAD, --iAlleleDepth  iEVA Allele Depth reported as Ref,Alt
      -iRR, --ReadRef       Number of reads mapping REF allele.
      -iRA, --ReadAlt       Number of reads mapping ALT allele.
      -iQR, --MeanRefQscore
                            Mean Q-score for REF allele.
      -iQA, --MeanAltQscore
                            Mean Q-score for ALT allele.
      -iRMQ, --RefMeanMappingQuality
                            Mean mapping quality for reads supporting REF allele.
      -iAMQ, --AltMeanMappingQuality
                            Mean mapping quality for reads supporting ALT allele.
      -iNCR, --NumberClippedReadsRef
                            Number of clipped reads mapping REF allele.
      -iNCA, --NumberClippedReadsAlt
                            Number of clipped reads mapping ALT allele.
      -iCR, --ClippedReadsRef
                            Fraction of clipped reads mapping REF allele.
      -iCA, --ClippedReadsAlt
                            Fraction of clipped reads mapping ALT allele.

## Required Arguments

- **-I, --input**
    - Path to input vcf file.

- **-Ref, --reference**
    - Path to reference file in fasta format.

- **-O, --outfile**
    - Path to output vcf file.

- **-L, --list**
    - **REQUIRED ONLY TO EXTRACT ATTRIBUTES FROM BAM FILE.**
    - Path to Bam list file, one file for each row, as explained in the [Usage](#how_to) section.
    - The sample name in the vcf *genotype* field must match the Read Group Sample Name `RG:SM` tag in the bam file header. A Bam Index file is required in the same path as the bam file.

## Reference Extraction Arguments

Attributes extracted from the `reference.fasta` file will be reported in the vcf *INFO* column.

### -SR, --SimpleRepeat

Enable *Simple Repeat Finder*. The `-SR` option checks if a variant in the *POS* field of the input vcf file falls into a repeated sequence for a specific [window size](#WS).

iEVA identifies two classes of repetitions:

1. **Simple Repeat Sequence**: In -SR, we defined 3 simple rules, **BY DEFAULT**, to identify a tandem simple repeated sequence:
   - *Number of repetitions* (N) >= 3. e.g., 3(CG) = CGCGCG.
   - *Nucleotides composing repeating unit* (dY) = 2. For example, 3(ATT) is a SR sequence, while 3(AATC) is not, due to 3 different nucleotides composing the repeating unit.
   - *Max length of repeating unit* (L) = 6. By default, the -SR option extracts only repeated sequences with L=6. For example, N(ATTAATT) is not considered a SR sequence.

   These default rules with the `-SR` option in iEVA are selected because these repeated sequences mostly influence base-calling and sequencing quality metrics [1]. However, for more complex repetitions (or implementing your own rules for repeated sequences), see the [SRList](#SRList) option.

2. **Homopolymer Sequence**: Defined as:
   - *Number of repetitions* (N) > 3.
   - *Number of nucleotides composing single repeated sequence unit* (dY) = 1. e.g., 4(A).

   For example, `GAGAGAGA` is a Repeated Sequence of a `GA` dinucleotide with N = 4. `AAAAAAA` is an Homopolymer of `A` nucleotide with N = 7.

> **Note:**
>
> The `-SR` argument has three output values: `SR=0` for None, `SR=1` for Repeated Sequence, and `SR=2` for Homopolymer Sequence.
>
> iEVA does not extract all repeated sequences in the window size but only verifies if a variant falls into it: only in this case, iEVA reports the repeated sequence class, otherwise SR=0.
>
> A complete list of Simple Repeat Sequences extracted by iEVA can be found in the iEVA package at `./iEVA/Supplementary/Simple-Repeated-Sequence-list.txt`.

### -SRList, --SimpleRepeatList

Path to the list containing simple repeated sequences (one per line) to extract with iEVA. Enabling the `-SRList` command, you ask iEVA to check if a variant falls in a repeated sequence reported in the list, neglecting the previously explained rules for simple repeated sequences in the [SimpleRepeat](#SimpleRepeat) option.

Users can choose to find any type of repeated sequence only if they provide a single repeating unit in the file. Compared to the [SimpleRepeat](#SimpleRepeat) option, you can define your own rules, using also sequences composed of more than two nucleotides (dY >= 2) and with a max length of repeating unit greater than 6 (L >= 6). For example, create a `RepSeq.txt` file as follows:

    CTT
    ACG
    ACCTG

In this case, iEVA finds all repeated sequences with the number of repetitions > 3 like *CTTCTTCTTCTTCTT*, *ACGACGACGACG* or *ACCTGACCTGACCTG*.

> **Warning:**
> 
> Using the `-SRList` command, iEVA will not extract the default repetitions reported in the `./iEVA/Supplementary/Simple-Repeated-Sequence-list.txt` file but only those given in the file list. However, iEVA always finds homopolymer sequences.
> If you want to extract also the default repetitions, you can add your repeating unit of interest in that file and use it to extract all requested repeated sequences.

### -SRL, --SimpleRepeatLength

Total length, as a number of nucleotides, of the repeated sequence in the `-SR` option in which the variant falls.

> **Note:**
> `GAGAGAGA`, for example, has `SRL=8`.

### -SRU, --SimpleRepeatUnit

Report the simple repeat unit composing the simple repeated sequence. For example, `ATTATTATT` has a simple repeat unit `SRU=ATT`.

### -PNC, --PseudoNucleotidesComposition

To achieve an accurate description of a sequence around a variant in the vcf input file, we implemented a *Pseudo nucleotide composition or PseKNC* for a specific [window size](#WS).

PseKNC allows to enclose nucleotide sequence composition using a number of attributes depending on kmer size.
In iEVA, we chose a dinucleotide profile (kmer size of 2), resulting in 16 attributes: *AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT*. Each dinucleotide combination is represented by a float number.

> **See Also:**
> Further details about Pseudo nucleotide composition can be found [here](#2).

> **Note:**
> `-PNC` is reported in the vcf as a list of 16 values: `PNC=AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT`.

### -RM, --RepeatMasker

Check if the variant falls in a lower case sequence in the `reference.fasta` file for a specific [window size](#WS).

Sequence in lower case is used for sequences identified by [RepeatMasker](http://www.repeatmasker.org/) as low-complexity or repetitive elements as reported in [Fasta format specification](https://www.ncbi.nlm.nih.gov/projects/SNP/snp_legend.cgi?legend=fasta).

> **Note:**
> `-RM` argument has two possible values: `RM=0` if the fasta sequence is upper case, `RM=1` if the fasta sequence is lower case.

### -GC, --gcContent

GC content of a sequence for a specific [window size](#WS).

### -VC, --VariantClass

Verify if the variant class is *SNV*, *Insertion*, or *Deletion*. Other variant classes (*e.g.,* InDel) are reported as *Sequence Alteration*.

> **Note:**
> `-VC` argument has 4 possible values: `VC=snv` for *SNV*, `VC=Ins` for *Insertion*, `VC=Del` for *Deletion*, and `VC=Alt` for *Sequence Alteration*.

### -WS, --WindowSize

Set the window size for [SimpleRepeat](#SimpleRepeat) and [SRList](#SRList), [gcContent](#gcContent), and [PseudoNucleotidesComposition](#PseudoNucleotidesComposition) arguments. *WS=[20-600]*, *Default=40*

## Bam Extraction Arguments

These arguments extract various information about reads mapping the sample variant position. Attributes are extracted from the flag field in the `.bam` file and will be annotated in the vcf genotype field for the bam reported in the [Required](#Required) option.

> **See Also:**
> Look at [bam specifications](https://samtools.github.io/hts-specs/SAMv1.pdf) for details about the bam format.

### -SBR, --StrandBiasReads

Fisher's exact test based on read orientation (R1+,R1-,R2+,R2-) to detect strand bias. Only for paired-end sequencing experiments.
This option excludes reads flagged as duplicates and unmapped.

### -UnMap, --UnMappedReads

Fraction of *unmapped reads* for a variant position.

### -MMQ, --MeanMappingQuality

Mean mapping quality for reads mapping the variant position.

### -MQ0, --MappingQualityZero

Fraction of reads mapping the variant position with Mapping Quality *MQ=0*.

### -NPA, --NotPrimaryAlignment

Fraction of reads mapping the variant position and flagged as *not primary alignment*.

### -SA, --SupplementaryAlignment

Fraction of reads mapping the variant position and flagged as *supplementary alignment*.

### -NP, --NotPairedReads

Fraction of reads mapping the variant position and flagged as *not paired*.

### -NPP, --NotProperPairedReads

Fraction of reads mapping the variant position and flagged as *not proper paired*.

### -AS, --AlignmentScore

Mean *Alignment Score* for reads mapping the variant position. Tested on bam obtained with [BWA](http://bio-bwa.sourceforge.net/).

### -TDR, --TotalDupReads

Fraction of total reads mapping the variant position and marked as *duplicate*.

## Genotype Extraction Arguments

Genotype extraction arguments allow for extracting information about the sample genotype analyzing `.bam` file. Most of the command line options in this section start with `-i` because they are computed by iEVA and are affected by `-SNVmbq`, `-SNVmmq`, `-INDELmbq`, and `-INDELmmq` arguments. They set SNV/InDels mapping quality and base quality threshold for all reads supporting the variant position. All reads below the threshold are filtered.

### -TDP, --TotalDPUnfilter

Total read depth for a variant position. No filter applied, including duplicate reads. Look at `-iDP` option for filtered read depth.

### -SNVmbq, --SNVMinBaseQuality

Minimum base quality threshold for base supporting *SNV* position. *SNVmbq=[0-66]*, *Default=12*

### -SNVmmq, --SNVMinMappingQuality

Minimum mapping quality threshold for reads supporting *SNV*. *SNVmmq=[0-60]*, *Default=30*

### -INDELmbq, --IndelMinBaseQuality

Minimum base quality threshold for base supporting *Insertions* and *Deletions*. *INDELmbq=[0-66]*, *Default=10*

### -INDELmmq, --IndelMinMappingQuality

Minimum mapping quality threshold for reads supporting *Insertions* and *Deletions*. *INDELmmq=[0-60]*, *Default=20*

### -iNDR, --NumberReadDupRef

Number of reads mapping the variant position and marked as duplicate on REF allele in the vcf *REF* column.

### -iNDA, --NumberReadDupAlt

Number of reads mapping the variant position and marked as duplicate on ALT allele in the vcf *ALT* column.

### -iDR, --DuplicateReference

Fraction of reads mapping the variant position and marked as duplicate on REF allele in the vcf *REF* column.

### -iDA, --DuplicateAlternate

Fraction of reads mapping the variant position and marked as duplicate on ALT allele in the vcf *ALT* column.

### -iDDup, --DeltaDuplicate

Difference between the fraction of duplicate reads for REF and ALT alleles (%REF - %ALT). A negative value shows preference in duplicate for REF allele, positive for ALT.

### -iDP, --iEvaDepth

iEva read depth for the variant position counted as reads supporting *REF* allele + *ALT* allele. Only reads flagged as *proper paired*, *proper mapped*, and *not duplicate* are used to extract `iDP`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

### -iAD, --iAlleleDepth

iEVA allele depth as the number of reads supporting *REF* and *ALT* alleles reported as `iAD=iRR,iRA`

### -iRR, --ReadRef

Number of reads supporting *REF* allele. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

### -iRA, --ReadAlt

Number of reads supporting *ALT* allele. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

### -iQR, --MeanRefQscore

Mean Q-score for base supporting *REF* allele. iEVA extracts `iQR` for reads supporting the REF allele depending on the `Variant Class` attribute as follows:

- **SNV**: Mean Q-score for all bases supporting REF allele.
- **Insertion**: Mean Q-score for all inserted bases supporting REF allele.
- **Deletion**: Mean Q-score for the base before and after deletion.

### -iQA, --MeanAltQscore

Mean Q-score for base supporting *ALT* allele. iEVA extracts `QA` for reads supporting the ALT allele depending on the `Variant Class` attribute as follows:

- **SNV**: Mean Q-score for all bases supporting ALT allele.
- **Insertion**: Mean Q-score for all inserted bases supporting ALT allele.
- **Deletion**: Mean Q-score for the base before and after deletion.

### -iRMQ, --RefMeanMappingQuality

Mean mapping quality for reads supporting *REF* allele. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

### -iAMQ, --AltMeanMappingQuality

Mean mapping quality for reads supporting *ALT* allele. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

### -iCR, --ClippedReadsRef

Fraction of clipped reads mapping *REF* allele in `iRR`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

### -iCA, --ClippedReadsAlt

Fraction of clipped reads mapping *ALT* allele in `iRA`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

### -iNCR, --NumberClippedReadsRef

Number of clipped reads mapping *REF* allele in `iRR`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

### -iNCA, --NumberClippedReadsAlt

Number of clipped reads mapping *ALT* allele in `iRA`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment*, or marked as *duplicate*.

# References

1. Kieleczawa J. Fundamentals of Sequencing of Difficult Templates—An Overview. Journal of Biomolecular Techniques : JBT. 2006;17(3):207-217. [PubMed PMID: 16870712](https://www.ncbi.nlm.nih.gov/pubmed/16870712).

2. Chen W, Lin H, Chou KC. Pseudo nucleotide composition or PseKNC: an effective formulation for analyzing genomic sequences. Mol Biosyst. 2015 Oct;11(10):2620-34. [PubMed PMID: 26099739](https://www.ncbi.nlm.nih.gov/pubmed/26099739).

# FAQ

### Could iEVA be considered a variant caller?

iEVA is not a variant caller because it does not include any genotyping algorithm. It simply needs as input a vcf file from variant calling and generates as output a new vcf file enriched with all extracted information.

### How does iEVA extract genotype information?

iEVA extracts genotype attributes like [iDP](#DP) or [iAD](#AD) by simply taking into account information about the variant position, *REF*, and *ALT* fields from the vcf input file. After that, iEVA iterates over all reads mapping the variant position in the bam file, looking for reads carrying the REF or ALT allele.

### Can I use any bam file or reference fasta file with iEVA?

Yes. You can use any bam or reference fasta file used in the variant calling analysis.

However, consider that some information can only be extracted if you made the respective analysis on the bam file. For example, you cannot extract information about the fraction of duplicate reads without a marking duplicate step on the bam file. In this case, not extracted options will be reported with a missing value `'.'` in the output vcf file.
In addition, bam analysis depends on pysam requirements.
