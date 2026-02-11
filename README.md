#  Kssdtree: an interactive Python package for phylogenetic analysis based on sketching technique
Kssdtree is a versatile Python package for phylogenetic analysis, offering three distinct pipelines: the Routine Pipeline, the Reference Subtraction Pipeline, and the GTDB-based Phylogenetic Placement Pipeline.

(1) Routine Pipeline: A general-purpose tool for phylogenetic analysis of user genomic data.
(2) Reference Subtraction Pipeline: Designed for intra-species phylogenomic analysis.
(3) GTDB-based Phylogenetic Placement Pipeline: Facilitates the search for similar genomes in the Genome Taxonomy Database (GTDB), conducting phylogenetic analysis alongside these genomes and positioning the input genomes within the entire prokaryotic tree of life.

Kssdtree also provides one-stop tree construction and visualization. It can handle DNA sequences in both fasta and fastq formats, whether gzipped or not. Additionally, Kssdtree is compatible with multiple platforms (Linux, MacOS, and Windows) and can be run using Jupyter notebooks.

# 1. Installation 
Kssdtree requires the Python 3 environment and the dependent packages pandas, pyqt5, ete3, and requests. If kssdtree is installed using the pip command, these dependencies will be installed automatically. For MacOS, it requires Python 3.8 or higher version. For Windows, it requires Python 3.6 version and the installation of the gzip tool (https://gnuwin32.sourceforge.net/packages/gzip.htm) for sequence decompression.
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

# 2. Command line Quick-Tutorial (Important!!!)
```
pip install kssdtree_cmd

kssdtree --help             
usage: kssdtree [-h] {routine,subtract,place} ...
subcommands:
  {routine,subtract,place}
    routine             Routine Pipeline
    subtract            Reference Subtraction Pipeline
    place               GTDB-based Placement Pipeline
```
## 2.1 Routine Pipeline

```
kssdtree routine --help
usage: kssdtree routine [-h] -i INPUT [-m METHOD] [-v VISUALIZE] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input folder (.fasta/.fastq files)
  -m METHOD, --method METHOD
                        method for constructing the tree, either 'nj' (NJ) or
                        'dnj' (DNJ)
  -v VISUALIZE, --visualize VISUALIZE
                        visualization mode, either 'r' (rectangle) or 'c'
                        (circle)
  -o OUTPUT, --output OUTPUT
                        output .newick file
                        
Example: kssdtree routine -i inputs -o output.newick
```
## 2.2 Reference Subtraction Pipeline

```
kssdtree subtract --help
usage: kssdtree subtract [-h] -i INPUT -r REFERENCE [-m METHOD] [-v VISUALIZE]
                         [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input folder (.fasta/.fastq files)
  -r REFERENCE, --reference REFERENCE
                        input a reference .fasta/.fastq files
  -m METHOD, --method METHOD
                        method for constructing the tree, either 'nj' (NJ) or
                        'dnj' (DNJ)
  -v VISUALIZE, --visualize VISUALIZE
                        visualization mode, either 'r' (rectangle) or 'c'
                        (circle)
  -o OUTPUT, --output OUTPUT
                        output .newick file
                        
Example: kssdtree subtract -i inputs -r reference.fasta -o output.newick
```
## 2.3 GTDB-based Phylogenetic Placement Pipeline

```
kssdtree place --help   
usage: kssdtree place [-h] -i INPUT [-m METHOD] [-v VISUALIZE] [-N NUMBER]
                      [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input .fasta/.fastq file or folder (.fasta/.fastq
                        files)
  -m METHOD, --method METHOD
                        method for constructing the tree, either 'nj' (NJ) or
                        'dnj' (DNJ)
  -v VISUALIZE, --visualize VISUALIZE
                        visualization mode, either 'r' (rectangle) or 'c'
                        (circle)
  -N NUMBER, --number NUMBER
                        maximum number of nearest reference genomes for
                        retrieving from GTDB database
  -o OUTPUT, --output OUTPUT
                        output folder

Example: kssdtree place -i test.fasta -o output
```

# 3. Python Quick-Tutorial 
## 3.1 Routine Pipeline

```
import kssdtree
kssdtree.quick(shuf_file='./shuf_files/L3K10.shuf', genome_files='your input genomes path', output='output.newick',  method='nj', mode='r')
```
## 3.2 Reference Subtraction Pipeline

```
import kssdtree
kssdtree.quick(shuf_file='./shuf_files/L3K10.shuf', genome_files='your input genomes path', output='output.newick', reference='your reference genome path', method='nj', mode='r')
```
## 3.3 GTDB-based Phylogenetic Placement Pipeline

```
import kssdtree
kssdtree.quick(shuf_file='./shuf_files/L3K9.shuf', genome_files='your input genomes path', output='your output path', database='gtdbr214', method='nj', mode='r', N=30)
```
For 'L3K10.shuf' and 'L3K9.shuf', if set parameter shuf_file='L3K10.shuf' or shuf_file='L3K9.shuf', kssdtree will download automatically them before performing quick or sketch function. If the automatic download fails, you can manually download them from https://zenodo.org/records/12699159 or current directory shuf_files. For other '*.shuf' files, such as 'L2K8.shuf', etc., kssdtree will be generated automatically by shuffle function. More usages about Kssdtree, please see Kssdtree documentation (https://kssdtree.readthedocs.io/en/latest).
# 4. How to cite
Hang Yang, Xiaoxin Lu, Jiaxing Chang, Qing Chang, Wen Zheng, Zehua Chen, Huiguang Yi, Kssdtree: an interactive Python package for phylogenetic analysis based on sketching technique, Bioinformatics, Volume 40, Issue 10, October 2024, btae566, https://doi.org/10.1093/bioinformatics/btae566