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

*iEVA* is a command line tool developed in python 2.7. To install iEVA digit in your terminal window: ::

    $ git clone https://github.com/Matteodigg/iEVA

or download it from https://github.com/Matteodigg/iEVA.

Enter in *iEVA* directory and execute: ::

    $ python setup.py install

Done! iEVA is installed in ``/usr/local/bin``. To use iEVA, simply type: ::

    $ iEVA -I path/to/input.vcf -O path/to/output.vcf -Ref path/to/reference.fasta -[Arg1] ...

Detailed explanation of iEVA command line arguments can be found in :ref:`iEVA arguments <options>` section.

.. seealso::

    If problem occurs installing python modules in :file:`REQUIREMENTS` file, please refers to `module homepage <http://pypi.python.org/pypi>`_ and try to install them separately using pypi::

     $ pip install module_name


Quickstart
----------
iEVA is a command line tool installed in ``usr/local/bin``. To use iEVA, simply digit on your command line: ::

	$ iEVA -I input.vcf -O output.vcf -L bam_list.txt -Ref reference.fasta -[Arg1] -[Arg2] ...

See Documentation for further information.

Documentation
-------------
A fully documentation about usage of iEVA tool can be found at http://ieva.readthedocs.io/en/latest/index.html


Status
------
Beta-version
