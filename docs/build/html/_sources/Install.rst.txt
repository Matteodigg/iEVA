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
