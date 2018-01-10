import os
from setuptools import setup

version = '0.1'

description = "A command-line tool for identification and Extraction of Variant Attributes."
cur_dir = os.path.dirname(__file__)
requirements = open(os.path.join(cur_dir, 'REQUIREMENTS')).read()

try:
    long_description = open(os.path.join(cur_dir, 'README.rst')).read()
except:
    long_description = description

setup(
    name = "iEVA",
    version = version,
    url = 'https://github.com/Matteodigg/iEVA',
    license = 'MIT',
    description = description,
    long_description = long_description,
    author = 'Matteo Di Giovannantonio',
    author_email = 'matteodeg@gmail.com',
    packages = ['iEVA'],
    install_requires = requirements,
 #    entry_points = {
	# 'console_scripts': ['iEVA = iEVA.iEVA:main']
	# },
     classifiers=[
         #'Development Status :: 3 - Alpha',
         #'Environment :: Console',
         'Programming Language :: Python',
     ],
)
