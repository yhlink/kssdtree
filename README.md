# Kssdree: an interactive Python package for phylogenetic analysis based on sketching technique

Kssdtree is a versatile Python package for phylogenetic analysis, including three different workflows: routine pipeline, reference subtraction pipeline and phylogenetic placement pipeline. The routine pipeline serves as a versatile tool for general-purpose phylogenetic analysis of users' genomic data.  The reference subtraction pipeline designs for intra-species phylogenetic analysis. The phylogenetic placement pipeline facilitates the search for similar genomes in the Genome Taxonomy Database (GTDB). It conducts phylogenetic analysis alongside these similar genomes and positions the input genomes within the entire prokaryotic tree of life. 

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
## 3.1 Routine pipeline
Demonstrate the routine pipeline of Kssdtree on ES29 dataset, including a multi-step implementation and a one-step implementation.
```
import kssdtree

# multi-step implementation
# step1、sketching 29 E.coli/Shigella strains (ES29) with k-mer substring space shuffled file (L3K10.shuf)
kssdtree.sketch(shuffle='L3K10.shuf', genomes='ES29', output='ES29_sketch')

# step2、calculating distance matrix 
kssdtree.dist(ref_sketch='ES29_sketch', qry_sketch='ES29_sketch', output='ES29.phylip', flag=0)

# step3、constructing tree with NJ
kssdtree.build(phylip='ES29.phylip', output='ES29.newick', method='nj')

# step4、visualizing tree 
kssdtree.visualize(newick='ES29.newick', mode='r')

# one-step implementation (recommendation)
kssdtree.quick(shuffle='L3K10.shuf', genomes='ES29', output='ES29.newick',  method='nj', mode='r')

```
![image](https://github.com/yhlink/kssdtree/blob/master/cases/case1.png)


## 3.2 Reference subtraction pipeline
Demonstrate the reference subtraction pipeline of Kssdtree on HG43 dataset, including a multi-step implementation and a one-step implementation.
```
import kssdtree

# multi-step implementation
# step1、sketching a human reference genome (hg38) and 43 human genomes (HG43) with L3K10.shuf
kssdtree.sketch(shuffle='L3K10.shuf', genomes='hg38.fa.gz', output='Ref_sketch', set_opt=True)
kssdtree.sketch(shuffle='L3K10.shuf', genomes='HG43', output='HG43_sketch', set_opt=True)

# step2、subtracting reference sketches (Ref_sketch) from input sketches (HG43_sketch)
kssdtree.subtract(ref_sketch='Ref_sketch', genomes_sketch=='HG43_sketch', output=='HG43_sub_sketch')

# step3、calculating distance matrix
kssdtree.dist(ref_sketch='HG43_sub_sketch', qry_sketch='HG43_sub_sketch', output='HG43.phylip', flag=0)

# step4、constructing tree with NJ
kssdtree.build(phylip='HG43.phylip', output='HG43.newick', method='nj')

# step5、visualizing tree 
kssdtree.visualize(newick='HG43.newick', taxonomy='HG43.txt', mode='r')
# Note: HG43.txt records names (accessions) and taxonomies of 43 human genomes.
HG01123    American Ancestry
HG01258    American Ancestry
HG01358    American Ancestry
...
# one-step implementation (recommendation)
kssdtree.quick(shuffle='L3K10.shuf', genomes='HG43', output='HG43.newick', reference='hg38.fa.gz', taxonomy='HG43.txt', method='nj', mode='r')

```
![image](https://github.com/yhlink/kssdtree/blob/master/cases/case2.png)

## 3.3 Phylogenetic placement pipeline
Demonstrate the phylogenetic placement pipeline of Kssdtree using an assembled prokaryotic genome, including a multi-step implementation and a one-step implementation.
```
import kssdtree

# multi-step implementation
# step1、sketching gtdbr214 (pre-sketching and hosting server) and an assembled prokaryotic genome - GCF_001228905.1 (PROK1) with L3K9.shuf
kssdtree.sketch(shuffle='L3K9.shuf', genomes='gtdbr214', output='gtdbr214_sketch', set_opt=True)
kssdtree.sketch(shuffle='L3K9.shuf', genomes='PROK1', output='PROK1_sketch', set_opt=True)

# step2、retrieving 30 closest sketches from GTDB (R214) sketches, combining PROK1_sketch
kssdtree.retrieve(ref_sketch='gtdbr214_sketch', qry_sketch='PROK1_sketch', output='PROK31_sketch', N=30)

# step3、calculating distance matrix
kssdtree.dist(ref_sketch='PROK31_sketch', qry_sketch='PROK31_sketch', output='PROK31.phylip', flag=0)

# step4、constructing tree with NJ
kssdtree.build(phylip='PROK31.phylip', output='PROK31.newick', method='nj')

# step5、visualizing tree 
kssdtree.visualize(newick='PROK31.newick', mode='r')

# one-step implementation (recommendation)
kssdtree.quick(shuffle='L3K9.shuf', genomes='PROK1', output='PROK31.newick', reference='gtdbr214', method='nj', mode='r', N=30)

```
![image](https://github.com/yhlink/kssdtree/blob/master/cases/case3.png)

More usages about Kssdtree, please see Kssdtree user manual (http://18.205.53.149:8000/kssdtree/kssdtree_user_manual.pdf).
 




