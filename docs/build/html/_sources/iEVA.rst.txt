.. _options:

iEVA arguments
==============

A full list of iEVA tool arguments with a short description can be accessed by typing in your command line :class:`-h` [- -help] option::

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



.. _Required:


Required arguments
------------------

* **-I, --input**

    Path to input vcf file

* **-Ref, --reference**

    Path to reference file in fasta format

* **-O, --outfile**

    Path to output vcf file

* **-L, --list**

    *REQUIRED ONLY TO EXTRACT ATTRIBUTES FROM BAM FILE.*

    |    Path to Bam list file, one file for each raw, as explained in :ref:`Usage<how_to>`.
    |    Sample name in vcf *genotype* field has to be the same of Read Group Sample Name :class:`RG:SM` tag in bam file header. Bam Index file required in the same path of bam file.



Reference extraction arguments
------------------------------

Attributes extracted from :file:`reference.fasta` file will be reported in vcf *INFO* column.

.. _SimpleRepeat:

-SR, --SimpleRepeat
^^^^^^^^^^^^^^^^^^^

Enable *Simple Repeat Finder*. :class:`-SR` option checks if a variant in *POS* field of input vcf file falls into a repeated sequence for a specific :ref:`window size <WS>`.

iEVA is able to identify two class of repetitions:

1. **Simple Repeat Sequence**: in -SR we defined 3 simple rules, **BY DEFAULT**, in order to identify a tandem simple repeated sequence:

* *Number of repetitions* (N) >= 3.  *e.g* 3(CG) = CGCGCG

* *Nucleotides composing repeating unit* (dY) = 2. For example, given a number of repetitions N=3, 3(ATT) is a SR sequence. 3(AATC), instead, is not a SR sequence because of 3 different nucleotides composing repeating unit.

* *Max length of repeating unit* (L) = 6. By default, -SR option extracts only repeated sequences with L=6. For example, N(ATTAATT) is not considered a SR sequence

We selected these rules by default with :class:`-SR` option in iEVA because these repeated sequences mostly influence base-calling and sequencing quality metrics [1]_. However, if you need to find more complex repetitions (or implementing your own rules for repeated sequences) look at :ref:`SRList` option.

2. **Homopolymer Sequence**: it is defined as:

* *Number of repetitions* (N) > 3

* *Number of nucleotides composing single repeated sequence unit* (dY) = 1.  *e.g* 4(A)

For example, `GAGAGAGA` is a Repeated Sequence of a `GA` dinucleotide with N = 4.
Instead, `AAAAAAA` is an Homopolymer of `A` nucleotide with N = 7.

.. note::
    :class:`-SR` argument has three output values: :class:`SR=0` for None, :class:`SR=1` for Repeated Sequence and :class:`SR=2` for Homopolymer Sequence.

    iEVA does not extract all repeated sequence in window size but only verify if a variant falls into it: only in this case, iEVA report repeated sequence class, otherwise SR=0.

    A complete list of Simple Repeat Sequence extracted by iEVA can be found in iEVA package at :file:`./iEVA/Supplementary/Simple-Repeated-Sequence-list.txt`


.. _SRList:


-SRList, --SimpleRepeatList
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Path to list containing simple repeated sequence (one for each raw) to extract with iEVA. Enabling :class:`-SRList` command, you ask to iEVA to check if a variant falls in a repeated sequence reported in list, neglecting previously explained rules for simple repeated sequences in :ref:`SimpleRepeat` option.

User can choose to find any type of repeated sequence only if provides single repeating unit in file. Compared to :ref:`SimpleRepeat` option, you can define your own rules, using also sequences composed by more than two nucleotides (dY >= 2) and with a max length of repeating unit greater than 6 (L >= 6). For example, try to create a :file:`RepSeq.txt` file as follow: ::

    CTT
    ACG
    ACCTG

In this case, iEVA finds all repeated sequences with number of repetitions > 3 like *CTTCTTCTTCTTCTT*, *ACGACGACGACG* or *ACCTGACCTGACCTG* 

.. warning::
    Using :class:`-SRList` command iEVA will not extract the default repetitions reported in :file:`./iEVA/Supplementary/Simple-Repeated-Sequence-list.txt` file but only those given in file list. However, iEVA always finds homopolymer sequences.
    If you want to extract also default repetitions, you can add your repeating unit of interest in that file and use it to extract all requested repeated sequences.


-SRL, --SimpleRepeatLength
^^^^^^^^^^^^^^^^^^^^^^^^^^

Total length, as number of nucleotides, of repeated sequence in :class:`SR` option in which variant falls.

.. note::
    `GAGAGAGA`, for example, has :class:`SRL=8`.


.. _PseudoNucleotidesComposition:


-SRU, --SimpleRepeatUnit
^^^^^^^^^^^^^^^^^^^^^^^^

Report simple repeat unit composing simple repeated sequence. For axample, `ATTATTATT` has a saimple repeat unit :class:`SRU=ATT`.

-PNC, --PseudoNucleotidesComposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to achieve an accurate description of a sequence around a variant in vcf input file, we implemented in iEVA a *Pseudo nucleotide composition or PseKNC* for a specific :ref:`window size <WS>`.

PseKNC allows to enclose nucleotide sequence composition using a number of attributes depending on kmer size.
In iEVA we choosed a dinucleotide profile (kmer size of 2), resulting in 16 attributes: *AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT*. Each dinucleotide combination is represented by a float number.

.. seealso::
    Further details about Pseudo nucleotide composition at [2]_.

.. note::
    :class:`-PNC` is reported in vcf as a list of 16 values: :class:`PNC=AA,AC,AG,AT,CA,CC,CG,CT,GA,GC,GG,GT,TA,TC,TG,TT`.


-RM, --RepeatMasker
^^^^^^^^^^^^^^^^^^^

Check if variant falls in a lower case sequence on :file:`reference.fasta` file for a specific :ref:`window size <WS>`.

Sequence in lower case is used for sequences identified by `RepeatMasker <http://www.repeatmasker.org/>`_ as low-complexity or repetitive elements as reported in `Fasta format specification <https://www.ncbi.nlm.nih.gov/projects/SNP/snp_legend.cgi?legend=fasta>`_

.. note::
    :class:`-RM` argument has two possible values. :class:`RM=0` if fasta sequence is upper case, :class:`RM=1` if fasta sequence is lower case.


.. _gcContent:


-GC, --gcContent
^^^^^^^^^^^^^^^^

GC content of a sequence for a specific `window size <WS>`_.


.. _VC:


-VC, --VariantClass
^^^^^^^^^^^^^^^^^^^

Verify if variant class is *SNV*, *Insertion* or *Deletion*. Other variant classes (*e.g* InDel) are reported as *Sequence Alteration*.

.. note::
    :class:`-VC` argument has 4 possible values. :class:`VC=snv` for *SNV*, :class:`VC=Ins` for *Insertion*, :class:`VC=Del` for *Deletion* and :class:`VC=Alt` for *Sequence Alteration*


.. _WS:


-WS, --WindowSize
^^^^^^^^^^^^^^^^^

Set window size for :ref:`-SR <SimpleRepeat>` and :ref:`-SRList <SRList>`, :ref:`-GC <gcContent>` and :ref:`-PNC <PseudoNucleotidesComposition>` arguments. *WS=[20-600]*, *Default=40*


.. _Genotype:


Bam extraction arguments
------------------------

These arguments extract several information about reads mapping sample variant position. Attributes are extracted from flag field in :file:`.bam` file and will be annotated in vcf genotype field for bam reported in :ref:`-L <Required>` option.

.. seealso::
    Look at `bam specifications <https://samtools.github.io/hts-specs/SAMv1.pdf>`_ for details about bam format.

-SBR, --StrandBiasReads
^^^^^^^^^^^^^^^^^^^^^^^

Fisher's exact test based on read orientation (R1+,R1-,R2+,R2-) to detect strand bias. Only for paired-end sequencing experiments.
This option excludes reads flagged as duplicates and unmapped.


-UnMap, --UnMappedReads
^^^^^^^^^^^^^^^^^^^^^^^

Fraction of *unmapped reads* for a variant position.


-MMQ, --MeanMappingQuality
^^^^^^^^^^^^^^^^^^^^^^^^^^

Mean mapping quality for reads mapping variant position.


-MQ0, --MappingQualityZero
^^^^^^^^^^^^^^^^^^^^^^^^^^

Fraction of reads mapping variant position with Mapping Quaility *MQ=0*.


-NPA, --NotPrimaryAlignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fraction of reads mapping variant position and flagged as *not primary alignment*.


-SA, --SupplementaryAlignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fraction of reads mapping variant position and flagged as *supplementary alignment*.


-NP, --NotPairedReads
^^^^^^^^^^^^^^^^^^^^^

Fraction of reads mapping variant position and flagged as *not paired*.


-NPP, --NotProperPairedReads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fraction of reads mapping variant position and flagged as *not proper paired*.


-AS, --AlignmentScore
^^^^^^^^^^^^^^^^^^^^^

Mean *Alignment Score* for reads mapping variant position. Tested on bam obtained with `BWA <http://bio-bwa.sourceforge.net/>`_


-TDR, --TotalDupReads
^^^^^^^^^^^^^^^^^^^^^

Fraction of total reads mapping variant position and marked as *duplicate*.


Genotype extraction arguments
-----------------------------

Genotype extraction arguments allow to extract information about sample genotype analyzing :file:`.bam` file.
Most of the command line options in this section start with ``-i`` because are computed by iEVA and are affected by :ref:`-SNVmbq <SNVmbq>`, :ref:`-SNVmmq <SNVmmq>`, :ref:`-INDELmbq <INDELmbq>` and :ref:`-INDELmmq <INDELmmq>` arguments. They set SNV/InDels mapping quality and base quality threshold for all reads supporting variant position. All reads below threshold are filtered.


-TDP, --TotalDPUnfilter
^^^^^^^^^^^^^^^^^^^^^^^

Total read depth for a variant position. No filter applied, including duplicate reads. Look at :class:`-iDP` option for filtered read depth.


.. _SNVmbq:


-SNVmbq, --SNVMinBaseQuality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Minimum base quality threshold for base supporting *SNV* position.

*SNVmbq=[0-66]*, *Default=12*


.. _SNVmmq:


-SNVmmq, --SNVMinMappingQuality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Minimum mapping quality threshold for reads supporting *SNV*.

*SNVmmq=[0-60]*, *Default=30*


.. _INDELmbq:


-INDELmbq, --IndelMinBaseQuality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Minimum base quality threshold for base supporting *Insertions* and *Deletions*.

*INDELmbq=[0-66]*, *Default=10*


.. _INDELmmq:


-INDELmmq, --IndelMinMappingQuality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Minimum mapping quality threshold for reads supporting *Insertions* and *Deletions*.

*INDELmmq=[0-60]*, *Default=20*


.. _NDR:


-iNDR, --NumberReadDupRef
^^^^^^^^^^^^^^^^^^^^^^^^^

Number of reads mapping variant position and marked as duplicate on REF allele in vcf *REF* column.


.. _NDA:


-iNDA, --NumberReadDupAlt
^^^^^^^^^^^^^^^^^^^^^^^^^

Number of reads mapping variant position and marked as duplicate on ALT allele in vcf *ALT* column.


.. _DR:


-iDR, --DuplicateReference
^^^^^^^^^^^^^^^^^^^^^^^^^^

Fraction of reads mapping variant position and marked as duplicate on REF allele in vcf *REF* column.


.. _DA:


-iDA, --DuplicateAlternate
^^^^^^^^^^^^^^^^^^^^^^^^^^

Fraction of reads mapping variant position and marked as duplicate on ALT allele in vcf *ALT* column.


.. _DDup:


-iDDup, --DeltaDuplicate
^^^^^^^^^^^^^^^^^^^^^^^^

Difference between fraction of duplicate reads for REF and ALT alleles (%REF - %ALT). A negative value shows preference in duplicate for REF allele, positive for ALT.


.. _DP:


-iDP, --iEvaDepth
^^^^^^^^^^^^^^^^^

iEva read depth for variant position counted as reads supporting *REF* allele + *ALT* allele. Only reads flagged as *proper paired*, *proper mapped* and *not duplicate* are used to extract :class:`iDP`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


.. _AD:


-iAD, --iAlleleDepth
^^^^^^^^^^^^^^^^^^^^

iEVA allele depth as number of reads supporting *REF* and *ALT* allele reported as :class:`iAD=iRR,iRA`


.. _RR:


-iRR, --ReadRef
^^^^^^^^^^^^^^^

Number of reads supporting *REF* allele. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


.. _RA:


-iRA, --ReadAlt
^^^^^^^^^^^^^^^

Number of reads supporting *ALT* allele. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


.. _QR:


-iQR, --MeanRefQscore
^^^^^^^^^^^^^^^^^^^^^

Mean Q-score for base supporting *REF* allele. iEVA extracts :class:`iQR` for reads supporting ref allele depending on :ref:`Variant Class <VC>` attribute as follow:

* **SNV**: Mean Q-score for all bases supporting REF allele

* **Insertion**: Mean Q-score for all inserted bases supporting REF allele

* **Deletion**: Mean Q-score for the base before and after deletion.


.. _QA:


-iQA, --MeanAltQscore
^^^^^^^^^^^^^^^^^^^^^

Mean Q-score for base supporting *ALT* allele. iEVA extracts :class:`QA` for reads supporting alt allele depending on :ref:`Variant Class <VC>` attribute as follow:

* **SNV**: Mean Q-score for all bases supporting ALT allele

* **Insertion**: Mean Q-score for all inserted bases supporting ALT allele

* **Deletion**: Mean Q-score for the base before and after deletion.


-iRMQ, --RefMeanMappingQuality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mean mapping quality for reads supporting *REF* allele. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


-iAMQ, --AltMeanMappingQuality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mean mapping quality for reads supporting *ALT* allele. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


.. _CR:


-iCR, --ClippedReadsRef
^^^^^^^^^^^^^^^^^^^^^^^

Fraction of clipped reads mapping *REF* allele in :ref:`-iRR <RR>`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


.. _CA:


-iCA, --ClippedReadsAlt
^^^^^^^^^^^^^^^^^^^^^^^

Fraction of clipped reads mapping *ALT* allele in :ref:`-iRA <RA>`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


.. _NCR:


-iNCR, --NumberClippedReadsRef
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Number of clipped reads mapping *REF* allele in :ref:`-iRR <RR>`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


.. _NCA:


-iNCA, --NumberClippedReadsAlt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Number of clipped reads mapping *ALT* allele in :ref:`-iRA <RA>`. iEVA filters all reads flagged as *unmapped*, *secondary alignment*, *supplementary alignment* or marked as *duplicate*.


References
----------

.. [1] Kieleczawa J. Fundamentals of Sequencing of Difficult Templates—An Overview. Journal of Biomolecular Techniques : JBT. 2006;17(3):207-217. `PubMed PMID: 16870712 <https://www.ncbi.nlm.nih.gov/pubmed/16870712>`_


.. [2] Chen W, Lin H, Chou KC. Pseudo nucleotide composition or PseKNC: an effective formulation for analyzing genomic sequences. Mol Biosyst. 2015 Oct;11(10):2620-34. `PubMed PMID: 19505943 <https://www.ncbi.nlm.nih.gov/pubmed/26099739>`_












