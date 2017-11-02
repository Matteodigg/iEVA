.. _ex:

Tutorial and examples
=====================

This is the tutorial for iEVA. Detailed :ref:`argument description <options>` is available.

Strating from a :file:`Test-norm.vcf` normalized file, a :file:`ucsc.hg19.fasta` reference and a :file:`bamlist.txt`, let's take a look at how to use iEVA.


Extracting attributes from fasta reference
------------------------------------------

Given an input vcf file from a variant caller containing these variants: ::

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

we verify if any variant falls in a repeated sequence, reporting its length and the single repeated unit for a window size of 300. In addition, we are going to check if variant falls in a lower case sequence for RepeatMasker: ::

    $ iEVA -v -I Test-norm.vcf -O iEVA-Test.vcf -Ref ucsc.hg19.fasta -SR -SRL -RM -SRU -WS 300

        Extracting attributes on: chr1

        Extracting attributes on: chr2

        Extracting attributes on: chr3

        Extracting attributes on: chr6

        Extracting attributes on: chr9

        Extracting attributes on: chr18

        Exatraction: Done

.. note::
    To extract only sequence features, argument :class:`-L` for bam list is not required.

Output file :file:`iEVA-Test.vcf` is now annotated as follow in *INFO* field: ::

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

You can observe that two variants have been found to be in a repeated sequence: a simple repeated sequence (:class:`SR=1` composed by TA nucleotides repetitions for a total lenght :class:`SRL=22`) and a homopolymer sequence (:class:`SR=2` of nucleotide T with total length :class:`SRL=16`), respectively.

Try to add other information about sequence nucleotide composition using *Pseudo Nucleotide Composition*, *GC* content and reporting *variant class*. Use as input the previous annotated output :file:`iEVA-Test.vcf` by simply adding requested arguments: ::

    $ iEVA -v -I Test-norm.vcf -O iEVA-Test.vcf -Ref ucsc.hg19.fasta -SR -SRL -RM -PNC -GC -WS 300 -VC

Alternatively, you can use as input the previous annotated output :file:`iEVA-Test.vcf` by simply adding requested arguments as follow: ::

    $ iEVA -v -I iEVA-Test.vcf -O iEVA-Test-PNC-GC-VC.vcf -Ref ucsc.hg19.fasta -PNC -GC -WS 300 -VC

        Extracting attributes on: chr1

        Extracting attributes on: chr2

        Extracting attributes on: chr3

        Extracting attributes on: chr6

        Extracting attributes on: chr9

        Extracting attributes on: chr18

        Exatraction: Done

In both cases, we have this output file: ::

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




Extracting Sample level informations
------------------------------------

As concern sample level extraction, we are going to test different attributes on genotype field ``Sample-01`` and ``Sample-02`` of vcf file.

Write a :file:`Bamlist.txt` file with bam path of those samples, one for each raw: ::

    path/to/Sample-01.bam
    path/to/Sample-02.bam

iEVA allows to check the overall quality of a called variants using information stored in sample bam file. For example, try to extract some reads information for each sample like fraction of unmapped reads, not paired reads, not proper paired reads, mapping quality 0 reads, strand bias in reads orientation and duplicate reads with :class:`-UnMap`, :class:`-NP`, :class:`-NPP`, :class:`-MQ0`, :class:`-SBR` and :class:`-TDR` options: ::

    $ iEVA -v -I iEVA-Test.vcf -O iEVA-Test-02.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -UnMap -NP -NPP -TDR -MQ0 -SBR

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr6

    Extracting attributes on: chr9

    Extracting attributes on: chr18

    Exatraction: Done


Output file has, in addition to sequence information previously extracted, requested sample specific attributes in *FORMAT* vcf field: ::

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

In the following example we are going to extract attributes on genotype field. For example, try to add iEVA read depth and iEVA allele depth with a base quality threshold of 20 (18 for InDels) and a mapping quality of 60 to extract the most informative reads. In addition, try to add mean *REF* q-score and mean *ALT* q-score: ::

    $ iEVA -v -I iEVA-Test-02.vcf -O iEVA-Test-03.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -iDP -iAD -iQR -iQA -SNVmbq 20 -INDELmbq 18 -SNVmmq 60 -INDELmmq 60

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr6

    Extracting attributes on: chr9

    Extracting attributes on: chr18

    Exatraction: Done

We get this output :file:`iEVA-Test-03.vcf` file: ::

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

For example, try to extract information about fraction of clipped reads for *REF* and *ALT* allele. Moreover, we want to check if duplicate reads supporting reference and alternate allele are well-balanced or bias affected. So, try to add :class:`-iCR` and :class:`-iCA` for clipped reads information and :class:`-iDR`, :class:`-iDA` and its difference :class:`-iDDup` for duplicate reads with default values for mapping and base quality threshold.

.. warning::
    Remember that optional arguments for :ref:`base quality threshold <SNVmbq>` and :ref:`mapping quality threshold <SNVmmq>` affect all the *allele-specific* arguments. In this case, default values will be used.

::

    $ iEVA -v -I iEVA-Test-03.vcf -O iEVA-Test-04.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -iCR -iCA -iDR -iDA -iDDup

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr6

    Extracting attributes on: chr9

    Extracting attributes on: chr18

    Exatraction: Done

This is the output :file:`iEVA-Test-04.vcf` vcf file: ::

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


Finally, try to extract attributes both from reference fasta file and bam file only for ``Sample-01``. In this case, simply remove ``Sample-02`` bam file from :file:`Bamlist.txt`. Missing values for ``Sample-02`` are indicated by a ``.`` for all requested optional arguments. ::

    $ iEVA -v -I Test-norm.vcf -O iEVA-Test-Only-Sample-01.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -SR -SRL -RM -GC -WS 300 -AS -UnMap -SA -DDup -iDP -iAD -iQR -iQA -iRMQ -iAMQ

    Sample 'Sample-02' will not be annotated. Missing bam file.

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr6

    Extracting attributes on: chr9

    Extracting attributes on: chr18

    Exatraction: Done

We get the following result: ::

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

A really important feature of iEVA is its flexibility and usage on more variant callers. For example, try to merge results coming out from different variant caller, obtaining a vcf file like this with different variant caller-specific annotations: ::

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

To retrieve sequence attributes and genotype information for both ``Sample-01`` and ``Sample-02`` use iEVA as showed in previous examples: ::

    $ iEVA -v -I merged-norm.vcf -O iEVA-merged.vcf -Ref ucsc.hg19.fasta -L Bamlist.txt -SR -SRL -RM -GC -WS 300 -AS -Unmap -SA -iDDup -MQ0 -iDP -iAD -iQR -iQA

    Extracting attributes on: chr1

    Extracting attributes on: chr2

    Extracting attributes on: chr3

    Extracting attributes on: chr9

    Exatraction: Done

As you can see, resulting file reports the input vcf file with iEVA attributes extracted for all variants, independently from the variant caller used: ::

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



