import kssd
import nj
import dnj
import toolutils
import os
import platform
import shutil
import time
import requests


def shuffle(k=None, s=None, l=None, o=None):
    kssd.write_dim_shuffle_file(k, s, l, o)


def sketch(shuf_file=None, genome_files=None, output=None, set_opt=None):
    if shuf_file is not None and genome_files is not None and output is not None:
        if not os.path.exists(genome_files):
            print('No such file or directory: ', genome_files)
            return False
        if set_opt is None:
            set_opt = False
        if not toolutils.allowed_file(genome_files):
            for filename in os.listdir(genome_files):
                if not toolutils.allowed_file(filename):
                    print('Genome format error for file:', filename)
                    return False
        if not os.path.exists(shuf_file):
            if shuf_file in ['L3K9.shuf', './L3K9.shuf', 'L3K10.shuf', './L3K10.shuf']:
                print('Downloading...', shuf_file)
                start_time = time.time()
                if shuf_file == 'L3K9.shuf' or shuf_file == './L3K9.shuf':
                    url = 'https://zenodo.org/records/12699159/files/L3K9.shuf?download=1'
                else:
                    url = 'https://zenodo.org/records/12699159/files/L3K10.shuf?download=1'
                headers = {'Accept-Encoding': 'gzip, deflate'}
                response = requests.get(url, headers=headers, stream=True)
                with open(shuf_file, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                end_time = time.time()
                if end_time - start_time > 200:
                    print("Network timeout, please manually download from https://zenodo.org/records/12699159")
                    return False
                print('Download finished: ', shuf_file)
            elif shuf_file in ['L2K8.shuf', 'L2K9.shuf', 'L3K11.shuf', './L2K8.shuf', './L2K9.shuf', './L3K11.shuf']:
                print('Shuffling...', shuf_file)
                file_name = shuf_file.split('.')[0]
                k = int(file_name[3:])
                if k == 11 or k == 10:
                    s = 6
                else:
                    s = 5
                l = int(file_name[1])
                shuffle(k, s, l, file_name)
                print('Shuffle finished: ', shuf_file)
            else:
                print('No such file or directory: ', shuf_file)
                return False
        print('Sketching...')
        start = time.time()
        if set_opt:
            kssd.dist_dispatch(shuf_file, genome_files, output, 1, 0, 0, '', '')
        else:
            kssd.dist_dispatch(shuf_file, genome_files, output, 0, 0, 0, '', '')
        end = time.time()
        print('Sketch spend time：%.2fs' % (end - start))
        print('Sketch finished!')
        return True
    else:
        print('Args error!!!')
        return False


def dist(genome_sketch=None, output=None, metric=None, flag=None):
    if genome_sketch is not None and output is not None:
        if not os.path.exists(genome_sketch):
            print('No such file or directory: ', genome_sketch)
            return False
        if flag is None:
            flag = 0
        if metric is None:
            metric = 'mash'

        print('Disting...')
        start = time.time()
        if '/' in output:
            output_dir = os.path.dirname(output)
            output_name = output.split('/')[-1]
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                print("Created directory:", output_dir)
        else:
            output_name = output
        if output_name.endswith(".phy") or output_name.endswith(".phylip"):
            if metric not in ['mash', 'aaf']:
                print('Metric type error, only supports mash or aaf distance')
                return False
            else:
                kssd.dist_dispatch(genome_sketch, output, genome_sketch, 2, 0, flag, metric, '')
                end = time.time()
                print('Dist spend time：%.2fs' % (end - start))
                print('Dist finished!')
                return True
        else:
            print('Output type error, only supports .phylip (.phy) format:', output_name)
            return False
    else:
        print('Args error!!!')
        return False


def combine(genome_sketch1=None, genome_sketch2=None, output=None):
    if genome_sketch1 is not None and genome_sketch2 is not None and output is not None:
        if not os.path.exists(genome_sketch1):
            print('No such file or directory: ', genome_sketch1)
            return False
        if not os.path.exists(genome_sketch2):
            print('No such file or directory: ', genome_sketch2)
            return False
        kssd.dist_dispatch(output, genome_sketch1, genome_sketch2, 3, 0, 0, '', '')
        return True


def getlist(genome_sketch=None, output=None):
    if genome_sketch is not None and output is not None:
        if not os.path.exists(genome_sketch):
            print('No such file or directory: ', genome_sketch)
            return False
        kssd.print_gnames(genome_sketch, output)
        return True


def retrieve(database=None, genome_sketch=None, output=None, N=None, method=None):
    if database is not None and genome_sketch is not None and output is not None:
        if method is None:
            method = 'nj'
        if method not in ['nj', 'dnj']:
            print('Only support nj and dnj methods!!!')
            return
        if not os.path.exists(genome_sketch):
            print('No such file or directory: ', genome_sketch)
            return False
        if database == 'gtdbr214':
            print('Retrieving...')
            start = time.time()
            if not os.path.exists(output):
                os.makedirs(output)
                print("Created directory:", output)
            else:
                print('Output path exist!!!')
                return False
            newick, accession_taxonomy = toolutils.upload_request(qry_sketch=genome_sketch, method=method, N=N)
            if newick is None:
                print('Server error!!!')
                return False
            with open(os.path.join(output, 'output.newick'), 'w') as f:
                f.write(newick)
            with open(os.path.join(output, 'output_accession_taxonomy.txt'), 'w') as f:
                for key, value in accession_taxonomy.items():
                    f.write("%s %s\n" % (key, value))
            end = time.time()
            print('Retrieve spend time：%.2fs' % (end - start))
            print('Retrieve finished!')
            return True
        else:
            print('Only support gtdbr214 database!!!')
            return False
    else:
        print('Args error!!!')
        return False


def build(phylip=None, output=None, method=None):
    if phylip is not None:
        if not os.path.exists(phylip):
            print('No such file or directory: ', phylip)
            return False
        if method is None:
            method = 'nj'
        if method not in ['nj', 'dnj']:
            print('Only support nj and dnj methods!!!')
            return False
        print('Building...')
        if '/' in output:
            output_dir = os.path.dirname(output)
            output_name = output.split('/')[-1]
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                print("Created directory:", output_dir)
        else:
            output_name = output
        if output_name.endswith(".nwk") or output_name.endswith(".newick"):
            start = time.time()
            if method == 'nj':
                state = nj.build(phylip, output)
            else:
                if platform.system() == 'Linux':
                    state = dnj.build(phylip, output, method)
                else:
                    state = nj.build(phylip, output)
            if state == 1:
                with open(output, 'r') as f:
                    lines = f.readlines()
                    newick = ''.join(lines)
                    newick = newick.replace('\n', '')
                with open(output, 'w') as f:
                    f.write(newick)
                end = time.time()
                print('Build spend time：%.2fs' % (end - start))
                print('Build finished!')
                return True
            else:
                print('phylip format error, Check that the phylip format is consistent with NJ or DNJ requirements!!!')
                return False
        else:
            print('Output type error, only supports .newick (.nwk) format:', output_name)
            return False
    else:
        print('Args error!!!')
        return False


def visualize(newick=None, taxonomy=None, mode=None):
    if newick is not None:
        if not os.path.exists(newick):
            print('No such file or directory: ', newick)
            return False
        if mode is None:
            mode = 'r'
        if taxonomy is not None and mode == 'c':
            print('Warning: this pipeline only support 'r' (rectangle) mode !!!')
            mode = 'r'
            toolutils.view_tree(newick, taxonomy, mode=mode)
        else:
            toolutils.view_tree(newick, taxonomy, mode=mode)
    else:
        print('Args error!!!')
        return False


def union(ref_sketch=None, output=None):
    if ref_sketch is not None and output is not None:
        if not os.path.exists(ref_sketch):
            print('No such file or directory: ', ref_sketch)
            return False
        kssd.sketch_union(ref_sketch, output)
        return True
    else:
        return False


def subtract(ref_sketch=None, genome_sketch=None, output=None, flag=None):
    if ref_sketch is not None and genome_sketch is not None and output is not None:
        if not os.path.exists(ref_sketch):
            print('No such file or directory: ', ref_sketch)
            return False
        if not os.path.exists(genome_sketch):
            print('No such file or directory: ', genome_sketch)
            return False
        if flag == 1:
            print('Subtracting...')
            start = time.time()
            kssd.sketch_operate(ref_sketch, output, genome_sketch)
            end = time.time()
            print('Subtract spend time：%.2fs' % (end - start))
            print('Subtract finished!')
            return True
        else:
            timeStamp = int(time.mktime(time.localtime(time.time())))
            print('Subtracting...')
            start = time.time()
            temp_txt = 'ref.txt'
            kssd.print_gnames(ref_sketch, temp_txt)
            nums = 0
            with open(temp_txt, 'r') as file:
                for line in file:
                    nums += 1
            if nums == 1:
                temp_union_sketch = ref_sketch
            else:
                temp_union_sketch = 'ref_union_sketch_' + str(timeStamp)
            r = union(ref_sketch=ref_sketch, output=temp_union_sketch)
            if not r:
                print('Union error!!!')
                return False
            kssd.sketch_operate(temp_union_sketch, output, genome_sketch)
            end = time.time()
            current_directory = os.getcwd()
            temp_dir = os.path.join(current_directory, temp_union_sketch)
            if platform.system() == 'Linux':
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
                if os.path.exists(temp_txt):
                    os.remove(temp_txt)
            else:
                pass
            print('Subtract spend time：%.2fs' % (end - start))
            print('Subtract finished!')
            return True
    else:
        print('Args error!!!')
        return False


def quick(shuf_file=None, genome_files=None, output=None, reference=None, database=None, method='nj', mode='r', N=0):
    if reference is None and database is None:
        if shuf_file is not None and genome_files is not None and output is not None:
            if toolutils.is_positive_integer(N) or toolutils.is_negative_integer(N):
                print("N must = 0 !!!")
                return False
            timeStamp = int(time.mktime(time.localtime(time.time())))
            temp_sketch = toolutils.rs() + '_sketch_' + str(timeStamp)
            temp_phy = toolutils.rs() + '_temp.phy'
            print('Step1...')
            if not toolutils.allowed_file(genome_files):
                num = toolutils.get_file_num(genome_files)
                if num == 1:
                    print('genome_files is a folder containing at least two . fasta or .fastq files!!!')
                    return False
            else:
                print('genome_files is a folder containing at least two . fasta or .fastq files, not a file!!!')
                return False
            s1 = sketch(shuf_file=shuf_file, genome_files=genome_files, output=temp_sketch, set_opt=False)
            if not s1:
                return False
            print('Step2...')
            if method == 'nj':
                s2 = dist(genome_sketch=temp_sketch, output=temp_phy, flag=0)
            else:
                s2 = dist(genome_sketch=temp_sketch, output=temp_phy, flag=1)
            if not s2:
                return False
            print('Step3...')
            s3 = build(phylip=temp_phy, output=output, method=method)
            if not s3:
                return False
            print('Step4...')
            visualize(newick=output, mode=mode)
            if platform.system() == 'Linux':
                current_directory = os.getcwd()
                temp_dir1 = os.path.join(current_directory, temp_sketch)
                temp_dir2 = os.path.join(current_directory, 'distout')
                if os.path.exists(temp_dir1):
                    shutil.rmtree(temp_dir1)
                if os.path.exists(temp_dir2):
                    shutil.rmtree(temp_dir2)
                if os.path.exists(temp_phy):
                    os.remove(temp_phy)
        else:
            print('Args error, please see https://kssdtree.readthedocs.io/en/latest!!!')
            return False
    elif reference is None and database == 'gtdbr214':
        if shuf_file is not None and genome_files is not None and output is not None:
            if not toolutils.is_positive_integer(N):
                print("N must > 0 !!!")
                return False
            if shuf_file != 'L3K9.shuf':
                print("shuf_file must be set to 'L3K9.shuf'")
                return False
            timeStamp = int(time.mktime(time.localtime(time.time())))
            qry_sketch = toolutils.rs() + '_sketch_' + str(timeStamp)
            s1 = sketch(shuf_file=shuf_file, genome_files=genome_files, output=qry_sketch, set_opt=True)
            if not s1:
                return False
            s2 = retrieve(database=database, genome_sketch=qry_sketch, output=output, N=N, method=method)
            if not s2:
                return False
            visualize(newick=os.path.join(output, 'output.newick'),
                      taxonomy=os.path.join(output, 'output_accession_taxonomy.txt'), mode=None)
            if platform.system() == 'Linux':
                current_directory = os.getcwd()
                temp_dir1 = os.path.join(current_directory, qry_sketch)
                temp_dir2 = os.path.join(current_directory, qry_sketch + '.zip')
                if os.path.exists(temp_dir1):
                    shutil.rmtree(temp_dir1)
                if os.path.exists(temp_dir2):
                    os.remove(temp_dir2)
        else:
            print('Args error, please see https://kssdtree.readthedocs.io/en/latest!!!')
            return False
    elif reference is None and database != 'gtdbr214':
        if shuf_file is not None and genome_files is not None and output is not None:
            if toolutils.is_positive_integer(N) or toolutils.is_negative_integer(N):
                print("N must = 0 !!!")
                return False
            if not os.path.exists(database):
                print('No such file or directory: ', database)
                return False
            if '/' in output:
                output_dir = os.path.dirname(output)
                output_name = output.split('/')[-1]
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                    print("Created directory:", output_dir)
            else:
                output_name = output
            if output_name.endswith(".newick"):
                timeStamp = int(time.mktime(time.localtime(time.time())))
                qry_sketch = toolutils.rs() + '_sketch_' + str(timeStamp)
                temp_combine_sketch = toolutils.rs() + '_combine_sketch_' + str(timeStamp)
                temp_phy = toolutils.rs() + '.phy'
                s1 = sketch(shuf_file=shuf_file, genome_files=genome_files, output=qry_sketch, set_opt=True)
                if not s1:
                    return False
                print('Step2...')
                combine(genome_sketch1=database, genome_sketch2=qry_sketch, output=temp_combine_sketch)
                if method == 'nj':
                    s2 = dist(genome_sketch=temp_combine_sketch, output=temp_phy, flag=0)
                else:
                    s2 = dist(genome_sketch=temp_combine_sketch, output=temp_phy, flag=1)
                if not s2:
                    return False
                print('Step3...')
                s3 = build(phylip=temp_phy, output=output, method=method)
                if not s3:
                    return False
                print('Step4...')
                getlist(genome_sketch=database, output='ref.txt')
                getlist(genome_sketch=qry_sketch, output='qry.txt')
                with open('ref.txt', 'r') as ref_file:
                    ref_lines = ref_file.readlines()
                with open('qry.txt', 'r') as qry_file:
                    qry_lines = qry_file.readlines()
                with open('ref_qry.txt', 'w') as result_file:
                    for line in ref_lines:
                        new_name = toolutils.rename_genome(line.strip())
                        result_file.write(new_name + '\tReference\n')
                    for line in qry_lines:
                        new_name = toolutils.rename_genome(line.strip())
                        result_file.write(new_name + '\tUnknown\n')
                os.remove('ref.txt')
                os.remove('qry.txt')
                os.remove(temp_phy)
                visualize(newick=output, taxonomy='ref_qry.txt', mode='r')
                if platform.system() == 'Linux':
                    current_directory = os.getcwd()
                    temp_dir1 = os.path.join(current_directory, qry_sketch)
                    if os.path.exists(temp_dir1):
                        shutil.rmtree(temp_dir1)
                    temp_dir2 = os.path.join(current_directory, temp_combine_sketch)
                    if os.path.exists(temp_dir2):
                        shutil.rmtree(temp_dir2)
            else:
                print('Output type error, only supports .newick format:', output_name)
                return False
        else:
            print('Args error, please see https://kssdtree.readthedocs.io/en/latest!!!')
            return False
    elif reference is not None and database is None:
        if shuf_file is not None and genome_files is not None and output is not None and method in ['nj', 'dnj']:
            if toolutils.is_positive_integer(N) or toolutils.is_negative_integer(N):
                print("N must = 0 !!!")
                return False
            if not toolutils.allowed_file(genome_files):
                num = toolutils.get_file_num(genome_files)
                if num == 1:
                    print('genome_files is a folder containing at least two . fasta or .fastq files!!!')
                    return False
            else:
                print('genome_files is a folder containing at least two . fasta or .fastq files, not a file!!!')
                return False
            timeStamp = int(time.mktime(time.localtime(time.time())))
            temp_reference_sketch = toolutils.rs() + '_ref_sketch_' + str(timeStamp)
            temp_genomes_sketch = toolutils.rs() + '_sketch_' + str(timeStamp)
            if not toolutils.allowed_file(reference):
                # cur_path = os.getcwd()
                # ref_path = os.path.join(cur_path, reference)
                num = toolutils.get_file_num(reference)
                if num == 1:
                    temp_union_sketch = temp_reference_sketch
                else:
                    temp_union_sketch = toolutils.rs() + '_ref_union_sketch_' + str(timeStamp)
            else:
                temp_union_sketch = temp_reference_sketch
            temp_subtract_sketch = toolutils.rs() + '_subtract_sketch_' + str(timeStamp)
            temp_phy = toolutils.rs() + '_temp.phy'
            print('Step1...')
            s1 = sketch(shuf_file=shuf_file, genome_files=reference, output=temp_reference_sketch, set_opt=True)
            if not s1:
                return False
            s2 = sketch(shuf_file=shuf_file, genome_files=genome_files, output=temp_genomes_sketch, set_opt=True)
            if not s2:
                return False
            print('Step2...')
            s3 = union(ref_sketch=temp_reference_sketch, output=temp_union_sketch)
            if not s3:
                return False
            s4 = subtract(ref_sketch=temp_union_sketch, genome_sketch=temp_genomes_sketch,
                          output=temp_subtract_sketch, flag=1)
            if not s4:
                return False
            print('Step3...')
            if method == 'nj':
                s5 = dist(genome_sketch=temp_subtract_sketch, output=temp_phy,
                          flag=0)
            else:
                s5 = dist(genome_sketch=temp_subtract_sketch, output=temp_phy,
                          flag=1)
            if not s5:
                return False
            print('Step4...')
            s6 = build(phylip=temp_phy, output=output, method=method)
            if not s6:
                return False
            print('Step5...')
            visualize(newick=output, mode=mode)
            if platform.system() == 'Linux':
                current_directory = os.getcwd()
                temp_dir1 = os.path.join(current_directory, temp_reference_sketch)
                temp_dir2 = os.path.join(current_directory, temp_genomes_sketch)
                temp_dir3 = os.path.join(current_directory, temp_union_sketch)
                temp_dir4 = os.path.join(current_directory, temp_subtract_sketch)
                temp_dir5 = os.path.join(current_directory, 'distout')
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
                if os.path.exists(temp_phy):
                    os.remove(temp_phy)
        else:
            print('Args error, please see https://kssdtree.readthedocs.io/en/latest!!!')
            return False

    else:
        print('Pipeline error, please see https://kssdtree.readthedocs.io/en/latest!!!')
        return False
