import random
import operator
import os
import requests
import json
import zipfile
import pandas as pd
import string

def allowed_file(filename):
    allowed_extensions = ['.fa', '.fa.gz', '.fasta', '.fasta.gz', '.fna', '.fna.gz', '.fastq', '.fastq.gz']
    return any(filename.endswith(ext) for ext in allowed_extensions)

def rs():
    letters = string.ascii_lowercase
    numbers = string.digits
    random_letters = ''.join(random.choice(letters) for i in range(6))
    random_numbers = ''.join(random.choice(numbers) for i in range(3))
    return random_letters + random_numbers

def is_positive_integer(num):
    if isinstance(num, int) and num > 0:
        return True
    else:
        return False

def is_negative_integer(num):
    if isinstance(num, int) and num < 0:
        return True
    else:
        return False

def randomcolor():
    colorArr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F']
    color = ""
    for i in range(6):
        color += colorArr[random.randint(0, 14)]
    return "#" + color


def str_insert(str_origin, pos, str_add):
    str_list = list(str_origin)
    str_list.insert(pos, str_add)
    str_out = ''.join(str_list)
    return str_out


def get_file_num(cwd):
    res = []
    get_dir = os.listdir(cwd)
    for i in get_dir:
        sub_dir = os.path.join(cwd, i)
        if os.path.isdir(sub_dir):
            get_file_num(sub_dir)
        else:
            ax = os.path.basename(sub_dir)
            res.append(ax)
    return len(res)


def decompress_zip(zip_path, dir_path):
    f = zipfile.ZipFile(zip_path, 'r')
    for file in f.namelist():
        f.extract(file, dir_path)
    f.close()


def upload_request(qry_sketch, method, N):
    zip_file = qry_sketch + '.zip'
    zip = zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED)
    for item in os.listdir(qry_sketch):
        zip.write(qry_sketch + os.sep + item)
    zip.close()
    url = "http://www.metakssdcoabundance.link/kssdtree/upload"
    header = {
        "kssdtree": 'upload'
    }
    data = {
        'method': method,
        "N": N
    }
    current_path = os.getcwd()
    files = {'file': open(os.path.join(current_path, zip_file), "rb")}
    res = requests.post(url=url, headers=header, data=data, files=files)
    response = res.text
    json_data = json.loads(response)
    if json_data['code'] == 200:
        newick = json_data['newick']
        accession_taxonomy = json_data['accession_taxonomy']
        return newick, accession_taxonomy
    else:
        return None, None

def view_tree(newick, taxonomy, mode):
    if taxonomy is None:
        from ete3 import PhyloTree, TreeStyle, NodeStyle, TextFace
        t = PhyloTree(newick, sp_naming_function=None)
        ts = TreeStyle()
        ts.mode = mode
        ts.show_leaf_name = False
        ts.show_scale = True
        ts.margin_bottom = 6
        ts.margin_top = 6
        ts.margin_left = 6
        ts.margin_right = 6
        ts.branch_vertical_margin = 10
        ts.extra_branch_line_type = 0
        ts.extra_branch_line_color = 'black'
        for node in t.traverse():
            nstyle = NodeStyle()
            if node.is_leaf():
                nstyle["fgcolor"] = "black"
                nstyle["shape"] = "circle"
            else:
                nstyle["fgcolor"] = "blue"
                nstyle["shape"] = "circle"
            node.img_style = nstyle
            branch_name_face = TextFace(node.dist, fsize=8, fgcolor='black', tight_text=False,
                                        bold=False)
            node.add_face(branch_name_face, column=0, position='branch-top')
        t.show(tree_style=ts)
    else:
        mode = 'r'
        from ete3 import PhyloTree, TreeStyle, NodeStyle, faces, AttrFace, CircleFace, TextFace
        def layout(node):
            if node.is_leaf():
                if node.species in species_colors:
                    C = CircleFace(radius=8, color=species_colors.get(node.species), style="circle")
                    C.opacity = 1
                    faces.add_face_to_node(C, node, 0, position="aligned")
                    N = AttrFace("name", fsize=14, fgcolor="black")
                    faces.add_face_to_node(N, node, 0)
                    style1 = NodeStyle()
                    style1["fgcolor"] = species_colors.get(node.species)
                    style1["size"] = 2
                    style1["vt_line_color"] = species_colors.get(node.species)
                    style1["hz_line_color"] = species_colors.get(node.species)
                    style1["vt_line_width"] = 1
                    style1["hz_line_width"] = 1
                    style1["vt_line_type"] = 0
                    style1["hz_line_type"] = 0
                    node.img_style = style1
                else:
                    N = AttrFace("name", fsize=20, fgcolor="red")
                    faces.add_face_to_node(N, node, 0)

        current_directory = os.getcwd()
        taxonomy_path = os.path.join(current_directory, taxonomy)
        if not os.path.exists(taxonomy_path):
            print('"The file taxonomy txt does not exist."')
            return

        all_accessions = []
        accession_taxonomy = {}
        with open(taxonomy_path, 'r') as file:
            for line in file:
                columns = line.split()
                column_1 = columns[0]
                column_2 = columns[1:]
                tempfile = ''
                for x in column_2:
                    tempfile = tempfile + x + ' '
                tempfile = tempfile[:-1]
                all_accessions.append(column_1)
                accession_taxonomy[column_1] = tempfile

        known_species = []
        for x in all_accessions:
            if accession_taxonomy[x] != 'Unknown':
                known_species.append(accession_taxonomy[x])

        temp_nwk = 'temp_kssdtree.newick'

        with open(newick, 'r') as f:
            lines = f.readlines()[0]
            for x in all_accessions:
                x_len = len(x)
                x_index = lines.index(x)
                loc_index = x_index + x_len + 8
                if x in accession_taxonomy.keys():
                    lines = str_insert(lines, loc_index, '[&&NHX:species=' + accession_taxonomy.get(x) + ']')
                else:
                    lines = str_insert(lines, loc_index, '[&&NHX:species=NewSpecies]')
        if os.path.exists(temp_nwk):
            os.remove(temp_nwk)
        with open(temp_nwk, 'w') as f:
            f.write(lines)
        unique_species = list(set(known_species))
        species_colors = {}

        temp_color1s = ['#9987ce', '#63b2ee']
        temp_color2s = ['#9987ce', '#63b2ee', '#76da91']
        for i in range(len(unique_species)):
            if len(unique_species) == 1:
                species_colors[unique_species[i]] = '#9987ce'
            elif len(unique_species) == 2:
                species_colors[unique_species[i]] = temp_color1s[i]
            elif len(unique_species) == 3:
                species_colors[unique_species[i]] = temp_color2s[i]
            else:
                species_colors[unique_species[i]] = randomcolor()
        species_colors = dict(sorted(species_colors.items(), key=operator.itemgetter(0)))
        t = PhyloTree(temp_nwk, sp_naming_function=None)
        for n in t.traverse():
            n.add_features(weight=random.randint(0, 50))
        ts = TreeStyle()
        ts.layout_fn = layout
        ts.mode = mode
        ts.show_leaf_name = False
        ts.show_branch_length = True
        ts.margin_bottom = 6
        ts.margin_top = 6
        ts.margin_left = 6
        ts.margin_right = 6
        ts.branch_vertical_margin = 10
        ts.extra_branch_line_type = 0
        ts.extra_branch_line_color = 'black'
        for species, color in species_colors.items():
            ts.legend.add_face(CircleFace(radius=8, color=color, style="circle"), column=0)
            ts.legend.add_face(TextFace(text=" " + species, fsize=14, fgcolor="black"), column=1)
        ts.legend_position = 4
        for node in t.traverse():
            if node.species == "NewSpecies":
                nst = NodeStyle()
                nst["bgcolor"] = "LightGrey"
                nst["fgcolor"] = "red"
                nst["shape"] = "circle"
                nst["vt_line_color"] = "red"
                nst["hz_line_color"] = "red"
                nst["vt_line_width"] = 2
                nst["hz_line_width"] = 2
                nst["vt_line_type"] = 0
                nst["hz_line_type"] = 0
                node.img_style = nst
                node.set_style(nst)
        # t.render("bubble_map.png", w=600, dpi=300, tree_style=ts)
        t.show(tree_style=ts)

def deal_gtdb_txt(temp_dist_output):
    data = pd.read_csv(temp_dist_output, delimiter='\t', header=None, skiprows=1)
    column_2 = data.iloc[:, 1]
    with open('new.txt', 'w') as file:
        for item in column_2:
            file.write(str(item) + '\n')
    with open('new.txt', 'r') as file:
        txt1_contents = file.read().splitlines()
    with open('gtdb.txt', 'r') as file:
        txt2_contents = file.read().splitlines()
    txt1_dict = {word: index for index, word in enumerate(txt1_contents)}
    result = []
    count = 0
    for word in txt2_contents:
        if word in txt1_dict:
            count += 1
            result.append(f"{count} {word}")
        else:
            result.append(f"0 {word}")
    with open('new_gtdb.txt', 'w') as file:
        for line in result:
            file.write(line + '\n')
    with open('new_gtdb.txt', 'r') as file:
        lines = file.readlines()
    result_lines = [line for line in lines if line.split()[0] != '0']
    with open('related_genomes_values.txt', 'w') as file:
        file.writelines(result_lines)


def deal_gtdb_phy(phy_filename):
    with open('related_genomes_values.txt', 'r') as file:
        lines = file.readlines()
    for i in range(len(lines)):
        columns = lines[i].split(' ')
        if len(columns) >= 2:
            file_name1 = columns[1].split('/')[-1][:15]
            columns[1] = 'GTDB_' + file_name1
            lines[i] = ' '.join(columns) + '\n'
    with open('modified_file.txt', 'w') as file:
        file.writelines(lines)
    new_accession_list = []
    new_gtdb_taxonomy_list = []
    with open('static/gtdbr214_accession_taxonomy.txt', 'r') as file:
        for line in file:
            accession, taxonomy = line.strip().split('\t')
            new_accession_list.append(accession)
            new_gtdb_taxonomy_list.append(taxonomy)
    accession_taxonomy = {}
    for i in range(len(new_accession_list)):
        accession_taxonomy[new_accession_list[i]] = new_gtdb_taxonomy_list[i]
    filename = 'modified_file.txt'
    new_filename = 'new_accession_taxonomy.txt'
    with open(filename, 'r') as file:
        with open(new_filename, 'w') as new_file:
            for line in file:
                columns = line.split()
                column_1 = columns[0]
                column_2 = columns[1]
                column_3 = accession_taxonomy.get(column_2)
                new_line = f"{column_1}\t{column_2}\t{column_3}\n"
                new_file.write(new_line)
    filename = 'new_accession_taxonomy.txt'
    known_species = []
    all_accessions = []
    dict1 = {}
    dict2 = {}
    with open(filename, 'r') as file:
        for line in file:
            columns = line.split()
            column_1 = columns[0]
            column_2 = columns[1]
            column_3 = columns[2:]
            tempfile = ''
            for x in column_3:
                tempfile = tempfile + x + ' '
            tempfile = tempfile[:-1]
            known_species.append(tempfile)
            dict1[column_1] = column_2
            dict2[column_2] = tempfile
    data = []
    with open(phy_filename, "r") as file:
        lines = file.readlines()
        for line in lines:
            data.append(line.strip().split())
    first_col = [row[0] for row in data[1:]]
    for i, item in enumerate(first_col):
        if item in dict1:
            first_col[i] = dict1[item]
        else:
            first_col[i] = first_col[i][:15]
    for i, item in enumerate(first_col):
        data[i + 1][0] = item
    for i in range(len(data)):
        if i == 0:
            pass
        else:
            all_accessions.append(data[i][0])
    if os.path.exists(phy_filename):
        os.remove(phy_filename)
    with open(phy_filename, "w") as file:
        for row in data:
            file.write(" ".join(row) + "\n")
