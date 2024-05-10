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
    sketch: sketching genomes into sketch and generating sketch files.
    :param shuffle: Kssdtree provide 'L3K9.shuf' and 'L3K10.shuf' files as input for genome sketching or decomposition. The default is 'L3K10.shuf'.
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
            if shuffle == 'L3K9.shuf':
                print('downloading...', shuffle)
                import http.client
                http.client.HTTPConnection._http_vsn = 10
                http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'
                url = 'http://18.205.53.149:8000/kssdtree/shuffle/' + shuffle
                headers = {
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.45 Safari/537.36'
                }
                r = requests.get(url, headers=headers, stream=True)
                with open(os.getcwd() + "\\" + shuffle, mode="wb") as f:
                    f.write(r.content)
                print('download finished!', shuffle)
            else:
                file_name = shuffle.split('.')[0]
                k = int(file_name[3:])
                if k == 10:
                    s = 6
                else:
                    s = 5
                l = int(file_name[1])
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


def dist(ref_sketch=None, qry_sketch=None, output=None, flag=None):
    """
    computing pairwise distances between reference and query genomes, and then generating a distance matrix in phylip format.
    :param ref_sketch: The folder path for sketch result files of reference genome files.
    :param qry_sketch: The folder path for sketch result files of query genome files.
    :param output: The output filename of distance matrix in phylip format.
    :param flag: 0 or 1. 0,1 is used to generate the distance matrix required by NJ (0 for diagonal elements) and DNJ (no diagonal elements) respectively.
    :return: null
    """
    if flag is None:
        flag = 0
    if ref_sketch is not None and qry_sketch is not None and output is not None:
        print('disting...')
        start = time.time()
        kssd.dist_dispatch(ref_sketch, output, qry_sketch, 2, 0, flag)
        end = time.time()
        print('dist spend time：%.2fs' % (end - start))
        print('dist finished!')
    else:
        print('args error!!!')


def retrieve(ref_sketch=None, qry_sketch=None, output=None, N=None):
    """
    retrieving N closest sketches from reference or GTDB (R214) sketches and combining query sketch files.
    :param ref_sketch: The folder path for sketch result files of reference genome files.
    :param qry_sketch: The folder path for sketch result files of query genome files.
    :param output: The output folder path for retrieve sketch result files of genome files.
    :param N: Max number of nearest reference genomes.
    :return: 0/1
    """
    if ref_sketch is not None and qry_sketch is not None and output is not None:
        if ref_sketch == 'gtdbr214_sketch':
            print('retrieving...')
            start = time.time()
            temp_related_sketch = 'related_sketch'
            reference = 'static/gtdbr214_sketch'
            kssd.dist_dispatch(reference, output, qry_sketch, 2, N, 3)
            kssd.print_gnames(reference, 'gtdb.txt')
            file_path1 = os.path.join(os.getcwd(), 'distout', 'distance.out')
            toolutils.deal_gtdb_txt(file_path1)
            kssd.grouping_genomes('new_gtdb.txt', reference, temp_related_sketch)
            kssd.dist_dispatch(output, qry_sketch, temp_related_sketch, 3, 0, 0)
            end = time.time()
            print('retrieve spend time：%.2fs' % (end - start))
            print('retrieve finished!')
        else:
            print('retrieving...')
            start = time.time()
            timeStamp = int(time.mktime(time.localtime(time.time())))
            temp_related_sketch = 'related_sketch_' + str(timeStamp)
            kssd.dist_dispatch(ref_sketch, output, qry_sketch, 2, N, 3)
            kssd.print_gnames(ref_sketch, 'gtdb.txt')
            file_path1 = os.path.join(os.getcwd(), 'distout', 'distance.out')
            toolutils.deal_gtdb_txt(file_path1)
            kssd.grouping_genomes('new_gtdb.txt', ref_sketch, temp_related_sketch)
            kssd.dist_dispatch(output, qry_sketch, temp_related_sketch, 3, 0, 0)
            end = time.time()
            file_path1 = 'new.txt'
            file_path2 = 'gtdb.txt'
            file_path3 = 'new_gtdb.txt'
            file_path4 = 'related_genomes_values.txt'
            file_path5 = 'modified_file.txt'
            file_path6 = 'new_accession_taxonomy.txt'
            if os.path.exists(file_path1):
                os.remove(file_path1)
            if os.path.exists(file_path2):
                os.remove(file_path2)
            if os.path.exists(file_path3):
                os.remove(file_path3)
            if os.path.exists(file_path4):
                os.remove(file_path4)
            if os.path.exists(file_path5):
                os.remove(file_path5)
            if os.path.exists(file_path6):
                os.remove(file_path6)
            print('retrieve spend time：%.2fs' % (end - start))
            print('retrieve finished!')
    else:
        print('args error!!!')


def build(phylip=None, output=None, method=None):
    """
    constructing tree with NJ or DNJ and generating tree in newick format.
    :param phylip: The distance matrix in phylip format.
    :param output: 'nj'(NJ) or 'dnj'(DNJ) method for constructing tree. The default is 'nj'.
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
    visualizing tree with ETE3 toolkit.
    :param newick: The tree in newick format.
    :param taxonomy:  The taxonomy information in txt format, which records the name (accession) of genome and its taxonomy. The default is None.
    :param mode: 'r'(rectangle) or 'c'(circle) mode for visualizing tree. The default is 'r'.
    :return: null
    """
    if mode is None:
        mode = 'r'
    if newick is not None:
        toolutils.view_tree(newick, taxonomy, mode=mode)
    else:
        print('args error!!!')


def union(ref_sketch=None, output=None):
    """
    :param sketch:
    :param output:
    :return:
    """
    if ref_sketch is not None and output is not None:
        kssd.sketch_union(sketch, output)
    else:
        print('args error!!!')


def subtract(ref_sketch=None, genomes_sketch=None, output=None, flag=0):
    """
    subtracting the ref_sketch from genomes_sketch and creating the remainder sketch files.
    :param ref_sketch: The folder path for reference sketch result files.
    :param genomes_sketch: The folder path for sketch result files of genome files.
    :param output: The output folder path for remainder sketch result files.
    :param flag: 0.
    :return: null
    """
    if ref_sketch is not None and genomes_sketch is not None and output is not None:
        if flag == 1:
            print('subtracting...')
            start = time.time()
            kssd.sketch_operate(ref_sketch, output, genomes_sketch)
            end = time.time()
            print('subtract spend time：%.2fs' % (end - start))
            print('subtract finished!')
        else:
            timeStamp = int(time.mktime(time.localtime(time.time())))
            temp_union_sketch = 'ref_union_sketch_' + str(timeStamp)
            print('subtracting...')
            start = time.time()
            union(ref_sketch=ref_sketch, output=temp_union_sketch)
            kssd.sketch_operate(temp_union_sketch, output, genomes_sketch)
            end = time.time()
            current_directory = os.getcwd()
            temp_dir = os.path.join(current_directory, temp_union_sketch)
            if platform.system() == 'Linux':
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
            else:
                pass
            print('subtract spend time：%.2fs' % (end - start))
            print('subtract finished!')
    else:
        print('args error!!!')


def quick(shuffle=None, genomes=None, output=None, reference=None, taxonomy=None, method='nj', mode='r', N=0):
    """
    simplifying pipeline and eliminating the necessity of organizing many intermediate files.
    :param shuffle: Kssdtree provide frequently-used 'L3K9.shuf' and 'L3K10.shuf' files as input for genome sketching or decomposition. The default is 'L3K10.shuf'. If you want to perform phylogenetic placement, you must use 'L3K9.shuf' file.
    :param genomes: The folder path for genome files. It supports the input of genome files in fasta/fastq formats.
    :param output: The output filename of tree in newick format.
    :param reference: The default is None, will perform the routine workflow. If you want to perform the reference subtraction workflow, you can set reference to the reference genome file or folder path. If you want to perform the phylogenetic placement, you must set reference to ‘gtdbr214’.
    :param taxonomy: The filename of taxonomy information in txt format, which records the name (accession) of genome and its taxonomy. The default is None.
    :param method: 'nj'(NJ) or 'dnj'(DNJ) method for constructing tree. The default is 'nj'.
    :param mode: 'r'(rectangle) or 'c'(circle) mode for visualizing tree. The default is 'r'.
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
                dist(ref_sketch=temp_sketch, qry_sketch=temp_sketch, output=temp_phy, flag=0)
            else:
                dist(ref_sketch=temp_sketch, qry_sketch=temp_sketch, output=temp_phy, flag=1)
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
    elif reference == 'gtdbr214' and taxonomy is None:
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
                union(ref_sketch=temp_reference_sketch, output=temp_union_sketch)
                subtract(ref_sketch=temp_union_sketch, genomes_sketch=temp_genomes_sketch,
                         output=temp_subtract_sketch, flag=1)
                print('step3...')
                if method == 'nj':
                    dist(ref_sketch=temp_subtract_sketch, qry_sketch=temp_subtract_sketch, output=temp_phy,
                         flag=0)
                else:
                    dist(ref_sketch=temp_subtract_sketch, qry_sketch=temp_subtract_sketch, output=temp_phy,
                         flag=1)
                print('step4...')
                build(phylip=temp_phy, output=output, method=method)
                if platform.system() == 'Linux':
                    pass
                else:
                    print('step5...')
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
