iEVA - a command line tool for Integration and Exatraction of Variant Attributes in NGS analysis
================================================================================================

.. _contact:

:Author: Matteo Di Giovannantonio and Mario Urtis
:Date: |today|
:Version: |version|
:Contact: matteodeg@gmail.com

iEVA is a python command line tool for integration and extraction of variant attributes. A lot of unused attributes are stored in *BAM* file and not reported from variant callers during variant calling analysis. iEVA expands the number of informative attributes obtained by parsing a *variant call format (.vcf)* file with a *reference.fasta* file and *.bam* file of sample stored in vcf.

iEVA works on python 2.7 using *pysam* and *pyfasta* modules. The current version has been tested with::

    pysam 0.9.1.4
    pyfasta 0.5.2

See the :ref:`Installation notes <installation>` for details.

Contents:
---------

.. toctree::
   :maxdepth: 2

   Intro
   Install
   Usage
   iEVA
   example
   FAQ




Indices and tables
==================

Contents:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


References
----------

:The python language:
      http://www.python.org


:The pysam homepage:
      http://pysam.readthedocs.io/en/latest/index.html


:The pyfasta homepage:
      http://pypi.python.org/pypi/pyfasta/

For any question or bug, please refer to :ref:`contact author <contact>` information.

