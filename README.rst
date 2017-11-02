iEVA
====

iEVA is a python command line tool for Integration and Extraction of Variant Attributes parsing a vcf file with a reference.fasta file and bam file for sample in input vcf file.


Installation
------------
To install iEVA, clone or download from git: ::

	$ git clone https://github.com/Matteodigg/iEVA

``cd`` in iEVA directory and execute: ::

	$ sudo python setup.py install


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

