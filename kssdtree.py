import kssd
import nj
import dnj

import toolutils
import os
import platform
import shutil
import time
import requests


def shuffle(k=8, s=5, l=2, o='default'):
    print('shuffling...')
    kssd.write_dim_shuffle_file(k, s, l, o)
    print('shuffle finished!')


def sketch(shuffle=None, genomes=None, output=None, set_opt=None):
    """
    reduces genomes into sketch and generates sketch files.
    :param shuffle: We provide ‘L3K9.shuf’ and ‘L3K10.shuf’ files as input for genome sketching or decomposition. The default is ‘L3K10.shuf’.
    :param genomes: The folder path for genome files. It supports the input of genome files in fasta/fastq formats.
    :param output: The output folder path for sketch result files of genome files.
    :param set_opt: Whether to do the set operation, default is False, if you want to do the set operation, you can set set_opt=True.
    :return: null
    """
    if set_opt is None:
        set_opt = False
    if shuffle is not None and genomes is not None and output is not None:
        current_directory = os.getcwd()
        shuf_file_path = os.path.join(current_directory, shuffle)
        if not os.path.exists(shuf_file_path):
            if shuffle == 'L3K9.shuf' or shuffle == 'L3K10.shuf':
                url = 'https://raw.githubusercontent.com/yhlink/kssdtree/master/shuffle_file/' + shuffle
                r = requests.get(url)
                with open(shuffle, 'wb') as f:
                    f.write(r.content)
            else:
                file_name = shuffle.split('.')[0]
                k = int(shuffle.split('.')[0][3:])
                if k == 10:
                    s = 6
                else:
                    s = 5
                l = int(shuffle[1])
                shuffle(k=k, s=s, l=l, o=file_name)
        print('sketching...')
        start = time.time()
        if set_opt:
            kssd.dist_dispatch(shuffle, genomes, output, 1, 0, 0)
        else:
            kssd.dist_dispatch(shuffle, genomes, output, 0, 0, 0)
        end = time.time()
        print('sketch spend time：%.2fs' % (end - start))
        print('sketch finished!')
    else:
        print('args error!!!')


def dist(ref_sketch=None, qry_sketch=None, output=None, N=None, flag=None):
    """
    computes pairwise distances between reference and query genomes, and then generates a distance matrix in phylip format.
    :param ref_sketch: The folder path for sketch result files of reference genome files.
    :param qry_sketch: The folder path for sketch result files of query genome files.
    :param output: The output filename of distance matrix in phylip format.
    :param N: Max number of nearest reference genomes. The default is 0 for computing pairwise distances between genomes.
    :param flag: 0, 1, or 2. 0,1 is used to generate the distance matrix required by nj (0 for diagonal elements) and dnj (no diagonal elements) respectively, and 2 does not generate the distance matrix.
    :return: null
    """
    if N is None:
        N = 0
    if flag is None:
        flag = 0
    if ref_sketch is not None and qry_sketch is not None and output is not None:
        print('disting...')
        start = time.time()
        kssd.dist_dispatch(ref_sketch, output, qry_sketch, 2, N, flag)
        end = time.time()
        print('dist spend time：%.2fs' % (end - start))
        print('dist finished!')
    else:
        print('args error!!!')


def build(phylip=None, output=None, method=None):
    """
    constructs tree with NJ or DNJ and generates tree in newick formatconstructs tree with NJ or DNJ and generates tree in newick format.
    :param phylip: The distance matrix in phylip format.
    :param output: ‘nj’ (NJ) or ‘dnj’ (DNJ) method for constructing tree. The default is ‘nj’.
    :param method: The output filename of tree in newick format.
    :return: null
    """
    if method is None:
        method = 'nj'
    if method not in ['nj', 'dnj']:
        print('method only support nj and dnj!!!')
        return
    if phylip is not None:
        print('building...')
        start = time.time()
        if output is None:
            output = 'kssdtree.newick'
        if method == 'nj':
            state = nj.build(phylip, output)
        else:
            if platform.system() == 'Linux':
                state = dnj.build(phylip, output, method)
            else:
                state = nj.build(phylip, output)
        if state == 1:
            nwk_path = os.path.join(os.getcwd(), output)
            with open(nwk_path, 'r') as f:
                lines = f.readlines()
                newick = ''.join(lines)
                newick = newick.replace('\n', '')
            with open(nwk_path, 'w') as f:
                f.write(newick)
            end = time.time()
            print('build spend time：%.2fs' % (end - start))
            print('build finished!')
    else:
        print('args error!!!')


def visualize(newick=None, taxonomy=None, mode=None):
    """
    visualizes tree with ETE3 toolkit.
    :param newick: The tree in newick format.
    :param taxonomy:  The taxonomy information in txt format, which records the name (accession) of genome and its taxonomy. The default is None.
    :param mode: ‘r’ (rectangle) or ‘c’(circle) mode for visualizing tree. The default is ‘r’.
    :return: null
    """
    if mode is None:
        mode = 'r'
    if newick is not None:
        toolutils.view_tree(newick, taxonomy, mode=mode)
    else:
        print('args error!!!')


def union(sketch=None, output=None):
    """
    duplicates those k-mers (integers) duplicated in sketch and creates union sketch files.
    :param sketch: The folder path for sketch result files of genome files.
    :param output: The output folder path for union sketch result files.
    :return: null
    """
    if sketch is not None and output is not None:
        print('unioning...')
        start = time.time()
        kssd.sketch_union(sketch, output)
        end = time.time()
        print('union spend time：%.2fs' % (end - start))
        print('union finished!')
    else:
        print('args error!!!')


def subtract(union_sketch=None, genomes_sketch=None, output=None):
    """
    subtracts the union_ sketch from genomes_sketch and creates the remainder sketch files.
    :param union_sketch: The folder path for union sketch result files.
    :param genomes_sketch: The folder path for sketch result files of genome files.
    :param output: The output folder path for remainder sketch result files.
    :return: null
    """
    if union_sketch is not None and genomes_sketch is not None and output is not None:
        print('subtracting...')
        start = time.time()
        kssd.sketch_operate(union_sketch, output, genomes_sketch)
        end = time.time()
        print('subtract spend time：%.2fs' % (end - start))
        print('subtract finished!')
    else:
        print('args error!!!')


def extract(sketch=None, output=None):
    """
    obtains genome filename from sketch and creates a txt file.
    :param sketch: The folder path for sketch result files of genome files.
    :param output: The output filename of genome information in txt format, which records the name of genome.
    :return: null
    """
    if sketch is not None and output:
        print('extracting...')
        start = time.time()
        kssd.print_gnames(sketch, output)
        end = time.time()
        print('extract spend time：%.2fs' % (end - start))
        print('extract finished!')
    else:
        print('args error!!!')


def group(sketch=None, txt=None, output=None):
    """
    groups part sketches from sketch and creates the grouped sketch files.
    :param sketch: The folder path for sketch result files of genome files.
    :param txt: The genome information in txt format, which records the name of genome.
    :param output: The output folder path for grouped sketch result files.
    :return: null
    """
    if sketch is not None and txt is not None and output is not None:
        print('grouping...')
        start = time.time()
        kssd.grouping_genomes(txt, sketch, output)
        end = time.time()
        print('group spend time：%.2fs' % (end - start))
        print('group finished!')
    else:
        print('args error!!!')


def combine(sketch1=None, sketch2=None, output=None):
    """
    combines two sketches and creates the combined sketch files.
    :param sketch1: The folder path for sketch result files of genome files.
    :param sketch2: The folder path for sketch result files of genome files.
    :param output: The output folder path for combined sketch result files.
    :return: null
    """
    if sketch1 is not None and sketch2 is not None and output is not None:
        print('combining...')
        start = time.time()
        kssd.dist_dispatch(output, sketch1, sketch2, 3, 0, 0)
        end = time.time()
        print('combine spend time：%.2fs' % (end - start))
        print('combine finished!')
    else:
        print('args error!!!')


def quick(shuffle=None, genomes=None, output=None, reference=None, taxonomy=None, method='nj', mode='r', N=0):
    """
    simplifies workflow and eliminates the necessity of organizing many intermediate files.
    :param shuffle: We provide frequently-used ‘L3K9.shuf’ and ‘L3K10.shuf’ files as input for genome sketching or decomposition. The default is ‘L3K10.shuf’. If you want to perform phylogenetic placement, you must use ‘L3K9.shuf’ file.
    :param genomes: The folder path for genome files. It supports the input of genome files in fasta/fastq formats.
    :param output: The output filename of tree in newick format.
    :param reference: The default is None, will perform the routine workflow. If you want to perform the reference subtraction workflow, you can set reference to the reference genome file or folder path. If you want to perform the phylogenetic placement, you must set reference to ‘gtdbr214’.
    :param taxonomy: The filename of taxonomy information in txt format, which records the name (accession) of genome and its taxonomy. The default is None.
    :param method: ‘nj’ (NJ) or ‘dnj’ (DNJ) method for constructing tree. The default is ‘nj’.
    :param mode: ‘r’ (rectangle) or ‘c’ (circle) mode for visualizing tree. The default is ‘r’.
    :param N: Max number of nearest reference genomes. The default is 0 for computing pairwise distances between genomes on routine and reference subtraction workflows. If you want to perform the phylogenetic placement, you can set N > 0.
    :return: null
    """
    if reference is None and taxonomy is None:
        if shuffle is not None and genomes is not None and output is not None:
            for filename in os.listdir(genomes):
                if not toolutils.allowed_file(filename):
                    print('Genome format error for file:', filename)
                    return 0
            timeStamp = int(time.mktime(time.localtime(time.time())))
            temp_sketch = genomes + '_sketch_' + str(timeStamp)
            temp_phy = 'temp.phy'
            print('step1...')
            sketch(shuffle=shuffle, genomes=genomes, output=temp_sketch, set_opt=False)
            print('step2...')
            if method == 'nj':
                dist(ref_sketch=temp_sketch, qry_sketch=temp_sketch, output=temp_phy, N=None, flag=0)
            else:
                dist(ref_sketch=temp_sketch, qry_sketch=temp_sketch, output=temp_phy, N=None, flag=1)
            print('step3...')
            build(phylip=temp_phy, output=output, method=method)
            if platform.system() == 'Linux':
                pass
            else:
                print('step4...')
                print('tree visualization finished!')
                visualize(newick=output, taxonomy=taxonomy, mode=mode)
            current_directory = os.getcwd()
            temp_dir1 = os.path.join(current_directory, temp_sketch)
            temp_dir2 = os.path.join(current_directory, 'distout')
            if platform.system() == 'Linux':
                if os.path.exists(temp_dir1):
                    shutil.rmtree(temp_dir1)
                if os.path.exists(temp_dir2):
                    shutil.rmtree(temp_dir2)
            else:
                pass
        else:
            print('args error!!!')
    elif reference == "gtdbr214" and taxonomy is None:
        if shuffle is not None and genomes is not None and output is not None and toolutils.is_positive_integer(N):
            if shuffle != 'L3K9.shuf':
                print("shuffle must be set to 'L3K9.shuf'")
                return 0
            for filename in os.listdir(genomes):
                if not toolutils.allowed_file(filename):
                    print('Genome format error for file:', filename)
                    return 0
            timeStamp = int(time.mktime(time.localtime(time.time())))
            temp_sketch = genomes + '_sketch_' + str(timeStamp)
            sketch(shuffle=shuffle, genomes=genomes, output=temp_sketch, set_opt=True)
            newick, accession_taxonomy = toolutils.upload_request(dir_name=temp_sketch, method=method, N=N)
            with open(output, 'w') as f:
                f.write(newick)
            with open('accession_taxonomy.txt', 'w') as f:
                for key, value in accession_taxonomy.items():
                    f.write("%s %s\n" % (key, value))
            if platform.system() == 'Linux':
                pass
            else:
                print('tree visualization finished!')
                visualize(newick=output, taxonomy='accession_taxonomy.txt', mode=None)
        else:
            print('args error or N<=0!!!')
    else:
        if shuffle is not None and genomes is not None and output is not None and method in ['nj', 'dnj']:
            timeStamp = int(time.mktime(time.localtime(time.time())))
            temp_reference_sketch = 'ref_sketch_' + str(timeStamp)
            temp_genomes_sketch = genomes + '_sketch_' + str(timeStamp)
            if not toolutils.allowed_file(reference):
                cur_path = os.getcwd()
                ref_path = os.path.join(cur_path, reference)
                num = toolutils.get_file_num(ref_path)
                if num == 1:
                    temp_union_sketch = temp_reference_sketch
                else:
                    temp_union_sketch = 'ref_union_sketch_' + str(timeStamp)
            else:
                temp_union_sketch = temp_reference_sketch
            temp_subtract_sketch = genomes + '_subtract_sketch_' + str(timeStamp)
            temp_phy = 'temp.phy'
            print('step1...')
            sketch(shuffle=shuffle, genomes=reference, output=temp_reference_sketch, set_opt=True)
            sketch(shuffle=shuffle, genomes=genomes, output=temp_genomes_sketch, set_opt=True)
            print('step2...')
            union(sketch=temp_reference_sketch, output=temp_union_sketch)
            print('step3...')
            subtract(union_sketch=temp_union_sketch, genomes_sketch=temp_genomes_sketch,
                     output=temp_subtract_sketch)
            print('step4...')
            if method == 'nj':
                dist(ref_sketch=temp_subtract_sketch, qry_sketch=temp_subtract_sketch, output=temp_phy, N=None, flag=0)
            else:
                dist(ref_sketch=temp_subtract_sketch, qry_sketch=temp_subtract_sketch, output=temp_phy, N=None, flag=1)
            print('step5...')
            build(phylip=temp_phy, output=output, method=method)
            if platform.system() == 'Linux':
                pass
            else:
                print('step6...')
                print('tree visualization finished!')
                visualize(newick=output, taxonomy=taxonomy, mode=mode)
            current_directory = os.getcwd()
            temp_dir1 = os.path.join(current_directory, temp_reference_sketch)
            temp_dir2 = os.path.join(current_directory, temp_genomes_sketch)
            temp_dir3 = os.path.join(current_directory, temp_union_sketch)
            temp_dir4 = os.path.join(current_directory, temp_subtract_sketch)
            temp_dir5 = os.path.join(current_directory, 'distout')
            if platform.system() == 'Linux':
                if os.path.exists(temp_dir1):
                    shutil.rmtree(temp_dir1)
                if os.path.exists(temp_dir2):
                    shutil.rmtree(temp_dir2)
                if os.path.exists(temp_dir3):
                    shutil.rmtree(temp_dir3)
                if os.path.exists(temp_dir4):
                    shutil.rmtree(temp_dir4)
                if os.path.exists(temp_dir5):
                    shutil.rmtree(temp_dir5)
            else:
                pass
        else:
            print('args error!!!')
