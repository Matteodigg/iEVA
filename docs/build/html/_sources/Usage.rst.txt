.. _how_to:

Usage
=====

iEVA requires three arguments in order to extract attributes from reference.fasta file::

    $ iEVA -I path/to/input.vcf -Ref path/to/reference.fasta -O path/to/out.vcf -[opts]

An additional argument is required to extract attributes from bam file:: 

    $ iEVA -I path/to/input.vcf -Ref path/to/reference.fasta -L path/to/Bam_list -O path/to/out.vcf -[opts]

:file:`Bam_list` file contains path to sample/s bam file in *--input* argument. One for each raw::

    path/to/Sample1.bam
    path/to/Sample2.bam
    .
    .
    .

Index :file:`.bai` is required to be in same path of :file:`.bam` file.

.. Note::

    If you are running for the first time :class:`pyfasta` module on :file:`reference.fasta`, it is going to take some time to write index file.


.. _Prereq:


Prerequisites
-------------

For a proper usage, iEVA meets the following prerequisites:

1. **Vcf normalization**:
    |    Actual iEVA release does not extract any allele-specific attribute on multi-allelic sites with multiple values in vcf *ALT* field.
    |    In addition, it is considered a good standard to *normalize* and *left-align* a vcf file after variant calling. Further details about *normalization* `here <http://genome.sph.umich.edu/wiki/Variant_Normalization>`_.

    You can normalize and split multi-allelic sites in vcf file, before using iEVA, with `bcftools norm <https://samtools.github.io/bcftools/bcftools.html>`_ using the following command line option: ::

        $ bcftools norm -m -both -f reference.fasta INPUT.vcf > OUTPUT.split-and-norm.vcf

2. **Bam header and sample name**:
    |    To use *bam extraction*, sample name in vcf *genotype* field need to be the same of :class:`RG:SM` tag in bam file.
    |    To modify bam header sample name, you can use `Picard tool <https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups>`_ or `Samtools <http://www.htslib.org/doc/samtools.html>`_.


.. warning::
    iEVA informs about extraction progress with :class:`--verbose` option. To use this argument, you have to **sort** vcf file.








