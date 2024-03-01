# Kssdree: an interactive Python package for phylogenetic analysis based on sketching technique

Kssdtree is a versatile Python package for phylogenetic analysis, including three different workflows: routine workflow, reference subtraction workflow and phylogenetic placement workflow. The routine workflow serves as a versatile tool for general-purpose phylogenetic analysis of users' genomic data.  The reference subtraction workflow designs for intra-species phylogenetic analysis. The phylogenetic placement workflow facilitates the search for similar genomes in the Genome Taxonomy Database (GTDB). It conducts phylogenetic analysis alongside these similar genomes and positions the input genomes within the entire prokaryotic tree of life. 

It also provides one-stop tree construction and visualization. It can handle DNA sequences of both fasta or fastq format, whether gzipped or not. Kssdtree can run on multiple platforms (Linux, Windows, and MacOS).

### Table of contents
1. [Requirements](#1-requirements)
2. [Installation](#2-installation)
3. [Quick Tutorial](#3-quick-tutorial)

# 1 Requirements
Linux, MacOS and Windows packages are built for Python3 and require the following dependent packages: pyqt5, ete3, requests. In the case Kssdtree is installed using pip command, these dependencies should be installed automatically. For Windows, it also requires the installation of the gzip tool for sequence decompression and gcc compiler.

# 2 Installation 
```

# For Windows
# create virtual environment
conda create --name=kssdtree python=3.x
# activate virtual environment
conda activate kssdtree
# install libpython and m2w64-toolchain 
conda install libpython m2w64-toolchain -c msys2
# install kssdtree
pip install kssdtree

#For Linux
# install kssdtree
pip install kssdtree

#For MacOS
# install gcc
brew install gcc
# install kssdtree (x86_64)
pip install kssdtree
# install kssdtree (arm64)
arch -x86_64 $(which python3) -m pip install kssdtree

```

# 3 Quick Tutorial
## 3.1 Routine workflow
```
import kssdtree

# The routine workflow on ES29 dataset.
kssdtree.quick(shuffle='L3K10.shuf', genomes='ES29', output='ES29.newick', 
reference=None, taxonomy=None, method='nj', mode='r', N=0)

```
--shuffle: We use the 'L3K10.shuf' file for k-mer substring space decomposition.

--genomes: The folder path of the ES29 dataset includes 29 E. coli/Shigella genomes.  

--output: We use 'ES29.newick' as the output filename of the tree in newick format.  

--reference: The default for the reference is None, which will perform the routine workflow.  

--taxonomy: The default for taxonomy is None, which will not provide the taxonomies of these genomes.  

--method: The default method is 'nj', which will use the NJ method to construct the tree.  

--mode: The default mode is 'r', which will use rectangle mode to visualize the tree.  

--N: The default value for N is 0, which will directly compute pairwise distances between genomes.  


## 3.2 Reference subtraction workflow
```
import kssdtree

# The reference subtraction workflow on HG43 dataset.
kssdtree.quick(shuffle='L3K10.shuf', genomes='HG43', output='HG43.newick', 
reference='hg38.fa.gz', taxonomy='HG43.txt', method='nj', mode='r', N=0)

```
--shuffle: We use the 'L3K10.shuf' file for k-mer substring space decomposition.  

--genomes: The folder path of the HG43 dataset includes 43 human genomes.  

--output: We use 'HG43.newick' as the output filename of the tree in newick format.  

--reference: We use 'hg38.fa.gz' (hg38 human genome) as a reference, which will perform the reference subtraction workflow.

--taxonomy: We include the taxonomy of these genomes in the HG43.txt file, which records the name (accession) of the genome and its taxonomy. If the taxonomy of these genomes is unknown, you can set the taxonomy to None.  
The contents of 'HG43.txt' are in the following format:  

HG01123     American Ancestry  

HG01258		American Ancestry  

HG01358		American Ancestry  

...  

--method: The default method is 'nj', which will use the NJ method to construct the tree.  

--mode: The default mode is 'r', which will use rectangle mode to visualize the tree.  

--N: The default value for N is 0, which will compute pairwise distances between genomes.  



## 3.3 Phylogenetic placement workflow
```
import kssdtree

# The phylogenetic placement workflow using an assembled genome GCF_001228905.1.
kssdtree.quick(shuffle='L3K9.shuf', genomes='BACT1', output='BACT11.newick', 
reference='gtdbr214', taxonomy=None, method='nj', mode='r', N=10)

```
--shuffle: We use the 'L3K9.shuf' file for k-mer substring space decomposition.  

--genomes: The folder path of the BACT1 dataset includes an assembled genome accession GCF_001228905.1.  

--output: We use 'BACT11.newick' as the output filename of the tree in newick format.  

--reference: We use 'gtdbr214' as a reference, which will perform the phylogenetic tree localization of prokaryotes workflow.  

--taxonomy: The default for taxonomy is None, which will not provide the taxonomies of these genomes.  

--method: The default method is 'nj', which will use the NJ method to construct the tree.  

--mode: The default mode is 'r', which will use rectangle mode to visualize the tree.  

--N: We set N=10 to query 10 nearest species from gtdbr214 reference.  




