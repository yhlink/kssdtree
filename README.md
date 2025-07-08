#  Kssdtree: an interactive Python package for phylogenetic analysis based on sketching technique
Kssdtree is a versatile Python package for phylogenetic analysis, offering three distinct pipelines: the Routine Pipeline, the Reference Subtraction Pipeline, and the GTDB-based Phylogenetic Placement Pipeline.

Routine Pipeline: A general-purpose tool for phylogenetic analysis of user genomic data.
Reference Subtraction Pipeline: Designed for intra-species phylogenomic analysis.
GTDB-based Phylogenetic Placement Pipeline: Facilitates the search for similar genomes in the Genome Taxonomy Database (GTDB), conducting phylogenetic analysis alongside these genomes and positioning the input genomes within the entire prokaryotic tree of life.
Kssdtree also provides one-stop tree construction and visualization. It can handle DNA sequences in both fasta and fastq formats, whether gzipped or not. Additionally, Kssdtree is compatible with multiple platforms (Linux, MacOS, and Windows) and can be run using Jupyter notebooks.
# 1. Installation 
Kssdtree requires the Python 3 environment and the dependent packages pandas, pyqt5, ete3, and requests. If Kssdtree is installed using the pip command, these dependencies will be installed automatically. For MacOS, it requires Python 3.8 or higher version. For Windows, it requires Python 3.6 version and the installation of the gzip tool(https://gnuwin32.sourceforge.net/packages/gzip.htm) for sequence decompression.
## 1.1 Linux

```
pip install kssdtree
```
## 1.2 MacOS

```
# (Optional) Install gcc (/opt/homebrew/bin/gcc-12) 
brew install gcc@12

# Create a virtual environment
conda create --name=kssdtree python=3.10

# Activate the virtual environment
conda activate kssdtree

# Install kssdtree
pip install kssdtree
```
## 1.3 Windows

```
# Create a virtual environment
conda create --name=kssdtree python=3.6.13

# Activate the virtual environment
conda activate kssdtree

# (Optional) Install libpython and m2w64-toolchain
conda install libpython m2w64-toolchain -c msys2

# Install kssdtree
pip install kssdtree
```
# 2. Quick-Tutorial 
## 2.1 Routine Pipeline

```
import kssdtree
kssdtree.quick(shuf_file='./shuf_files/L3K10.shuf', genome_files='your input genomes path', output='output.newick',  method='nj', mode='r')
```
## 2.2 Reference Subtraction Pipeline

```
import kssdtree
kssdtree.quick(shuf_file='./shuf_files/L3K10.shuf', genome_files='your input genomes path', output='output.newick', reference='your reference genome path', method='nj', mode='r')
```
## 2.3 GTDB-based Phylogenetic Placement Pipeline

```
import kssdtree
kssdtree.quick(shuf_file='./shuf_files/L3K9.shuf', genome_files='your input genomes path', output='your output path', database='gtdbr214', method='nj', mode='r', N=30)
```
More usages about Kssdtree, please see Kssdtree documentation (https://kssdtree.readthedocs.io/en/latest).
# 3. How to cite
Hang Yang, Xiaoxin Lu, Jiaxing Chang, Qing Chang, Wen Zheng, Zehua Chen, Huiguang Yi, Kssdtree: an interactive Python package for phylogenetic analysis based on sketching technique, Bioinformatics, Volume 40, Issue 10, October 2024, btae566, https://doi.org/10.1093/bioinformatics/btae566