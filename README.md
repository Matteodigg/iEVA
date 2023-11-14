# iEVA

iEVA is a python command line tool for the Integration and Extraction of Variant Attributes.

## Problem

It is a good practice in DNAseq bioinformatics pipelines to include multiple variant callers to reach maximum sensitivity. However, merging information from all these variant callers in a harmonised and clear way is quite tricky.

Existing tools such as bcftools or GATK CombineVariants can merge multiple vcf files but only the annotation from one of them will be reported in the final merged file. Moreover, bioinformaticians and data scientists might be interested in getting additional information for each variant in the vcf files from different sources. For example, one might be interested in getting the GC content around the variant or the sequence conservation score. This information can only be obtained by parsing different datasets and the reference fasta file with the input vcf files and, eventually, the BAM file of the sample used to generate the vcf.

## Solution

iEVA is a python based command line tool that enables users to:

1. Merge multiple VCF files from different variant callers on both single- and multi-sample level.
2. Include external sequence-specific information including features from fasta reference file such as conservation score of the sequence, GC content, homopolymer content etc.
3. Generate a merged VCF file containing all the annotation from EACH input vcf file, conserving all the useful annotation derived from each variant caller, in addition to additional information explained above.

.. _installation:

Installation
============

*iEVA* is a command line tool developed in python 2.7. To install iEVA digit in your terminal window:

    $ git clone https://github.com/Matteodigg/iEVA

or download it from https://github.com/Matteodigg/iEVA.

Enter in *iEVA* directory and execute:

    $ python setup.py install

Done! iEVA is installed in ``/usr/local/bin``. To use iEVA, simply type:

    $ iEVA -I path/to/input.vcf -O path/to/output.vcf -Ref path/to/reference.fasta -[Arg1] ...

Detailed explanation of iEVA command line arguments can be found in :ref:`iEVA arguments <options>` section.

.. seealso:

    If problem occurs installing python modules in :file:`REQUIREMENTS` file, please refers to `module homepage <http://pypi.python.org/pypi>`_ and try to install them separately using pypi:

     $ pip install module_name

Introduction
============

iEVA is a python command line tool that expands the number of informative attributes that are not commonly provided by variant callers. It needs as input a ``.vcf`` file from variant calling and generates as output a new vcf file enriched with all extracted information.

iEVA attributes are annotated in *INFO* and *FORMAT* field of vcf output file (see https://samtools.github.io/hts-specs/VCFv4.2.pdf for *vcf* specification). iEVA extracts information for variant in input vcf file from two different sources:

1. **Reference.fasta file**. Nucleotide sequence composition around a called variant is extracted from the fasta reference file. The annotations are reported in vcf *INFO* field.

2. **BAM file**. Information about reads mapping variant position, for sample in vcf input file, is extracted from bam file. The annotations are reported in vcf *FORMAT* field.

iEVA enriches vcf file with several attributes. With iEVA you can:

* Report the reference sequence quality in terms of G-C content, repeated sequences, pseudo nucleotide composition and more over. 
* Report genotype information for samples in vcf file. You can extract, for each sample, useful information like fraction of unmapped reads, not paired reads, duplicate reads, etc.\.\.
* Report *allele-specific* information, by evaluating, for example, number of reads supporting reference or alternate allele given a specific mapping quality or base quality threshold. These attributes help to identify any bias affecting called variants.

Main advatages of iEVA are:

* Output file in a vcf format file. iEVA simply adds the extracted attributes, maintaining original tags in input vcf file.
* Extraction is independent of the variant caller used in variant calling analysis. To a proper usage of iEVA look at :ref:`Usage <how_to>` and :ref:`Prerequisites <Prereq>`.
* iEVA is not a variant caller, for this reason is a safe-time approach. However, its computational time increases depending both on the number of requested sample for extraction, target length and mean target coverage.

See :ref:`iEVA arguments <options>` for details about iEVA options and :ref:`Tutorial and examples <ex>` section.

One of the strengths of iEVA consists in its simplicity of integration in bioinformatics pipeline, as showed in this figure:

[workflow](/docs/source/Workflow-iEVA.PNG)

Quickstart
----------
iEVA is a command line tool installed in ``usr/local/bin``. To use iEVA, simply digit on your command line:

	$ iEVA -I input.vcf -O output.vcf -L bam_list.txt -Ref reference.fasta -[Arg1] -[Arg2] ...

See Documentation for further information.

Documentation
-------------
A fully documentation about usage of iEVA tool can be found at http://ieva.readthedocs.io/en/latest/index.html


Status
------
Beta-version
