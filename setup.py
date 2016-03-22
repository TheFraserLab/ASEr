"""ASEr Setup Script."""
from setuptools import setup, find_packages
from codecs import open
from os import path
from os import listdir

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Generate a list of python scripts
scpts = []
for i in listdir(here + '/bin'):
    if i.endswith('.py'):
        scpts.append('bin/' + i)

setup(
    name='ASEr',
    version='0.1',
    description='Get ASE counts from BAMs or raw fastq data -- repackage of pipeline by Carlo Artieri ',
    long_description=long_description,
    url='https://github.com/MikeDacre/ASEr',
    author='Michael Dacre',
    author_email='mike.dacre@gmail.com',
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Beta',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'Operating System :: Linux',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],

    keywords='ASE allele-specific expression RNA-seq fastq BAM SAM SNP',

    install_requires=['pybedtools', 'pysam'],
    scripts=scpts,
    packages=['ASEr']

)
