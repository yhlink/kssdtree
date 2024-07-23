#include <Python.h>
#include "kssdheaders/command_dist.h"
#include "kssdheaders/command_dist_wrapper.h"
#include "kssdheaders/command_shuffle.h"
#include "kssdheaders/global_basic.h"
#include "kssdheaders/command_set.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

#include <time.h>

#ifdef _WIN32

#include <malloc.h>

#elif __linux__
#include <malloc.h>
#elif __APPLE__
#include <stdlib.h>
#endif

dim_shuffle_stat_t dim_shuffle_stat = {
        0,
        8,
        5,
        2,
};
char shuf_out_file_prefix[PATHLEN] = "./default";

char *rename_file(char *filename, char suffix[]) {
    char *p1 = strrchr(filename, '/');
    p1 += 1;
    int len = strlen(p1);
    int suffix_len = strlen(suffix);
    char *newfile = (char *) malloc(PATHLEN * sizeof(char *));
    newfile = strncpy(newfile, p1, len - suffix_len);
    newfile[len - suffix_len] = '\0';
    return newfile;
}


int endsWith(char *str, char *suffix) {
    int str_len = strlen(str);
    int suffix_len = strlen(suffix);
    if (str_len < suffix_len) {
        return 0;
    }
    return (strcmp(str + str_len - suffix_len, suffix) == 0);
}


int endsWithAny(char *str, char *suffixes[], int num_suffixes) {
    for (int i = 0; i < num_suffixes; i++) {
        if (endsWith(str, suffixes[i])) {
            return 1;
        }
    }
    return 0;
}

int create_matrix(char *input_name, char *output_name, int flag) {
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();
    FILE *fp = fopen(input_name, "r");
    if (fp == NULL) {
        fprintf(stderr, "can't open file %s", input_name);
        return 0;
    }
    FILE *fo = fopen(output_name, "w");
    char line[PATHLEN];
    char temp[PATHLEN];
    char seq1[PATHLEN];
    char seq2[PATHLEN];
    double distance;
    int num_seqs = 0;
    fgets(line, PATHLEN, fp);
    char pre_seq[PATHLEN];
    fgets(line, PATHLEN, fp);
    sscanf(line, "%s", temp);
    strcpy(pre_seq, temp);
    rewind(fp);
    fgets(line, PATHLEN, fp);
    char **seq_names = NULL;
    int max_num = 1000;
    seq_names = (char **) malloc(max_num * sizeof(char *));
    if (seq_names == NULL) {
        fprintf(stderr, "memory malloc error!\n");
    }
    while (fgets(line, PATHLEN, fp) != NULL) {
        sscanf(line, "%s %s", seq1, seq2);
        if (strcmp(pre_seq, seq1) != 0) {
            break;
        }
        if (num_seqs >= max_num) {
            max_num *= 2;
            seq_names = (char **) realloc(seq_names, max_num * sizeof(char *));
            if (seq_names == NULL) {
                fprintf(stderr, "memory malloc error!\n");
            }
        }
        seq_names[num_seqs] = (char *) malloc((strlen(seq2) + 1) * sizeof(char));
        if (seq_names[num_seqs] == NULL) {
            fprintf(stderr, "memory malloc error!\n");
        }
        strcpy(seq_names[num_seqs], seq2);
        num_seqs++;
    }
    double **distances = malloc(num_seqs * sizeof(double *));
    if (flag == 0) {
        for (int i = 0; i < num_seqs; i++) {
            distances[i] = malloc((i + 1) * sizeof(double));
            for (int j = 0; j <= i; j++)
                distances[i][j] = 0.0;
        }
    } else {
        for (int i = 0; i < num_seqs; i++) {
            distances[i] = malloc(i * sizeof(double));
            for (int j = 0; j < i; j++)
                distances[i][j] = 0.0;
        }
    }
    rewind(fp);
    fgets(line, PATHLEN, fp);
    int i = 0;
    int j = 0;
    if (flag == 0) {
        while (fgets(line, PATHLEN, fp) != NULL) {
            sscanf(line, "%*s %*s %*s %*s %lf", &distance);
            if (j <= i) {
                distances[i][j] = distance;
            }
            i += 1;
            if (i == num_seqs && j < num_seqs) {
                i = 0;
                j += 1;
            }
        }
    } else {
        while (fgets(line, PATHLEN, fp) != NULL) {
            sscanf(line, "%*s %*s %*s %*s %lf", &distance);
            if (j < i) {
                distances[i][j] = distance;
            }
            i += 1;
            if (i == num_seqs && j < num_seqs) {
                i = 0;
                j += 1;
            }
        }
    }

    for (int i = 0; i < num_seqs; i++) {
        char *suffixes[] = {".fasta.gz", ".fasta", ".fastq.gz", ".fastq", ".fna.gz", ".fna", ".fa.gz", ".fa"};
        int num_suffixes = 8;
        if (endsWithAny(seq_names[i], suffixes, num_suffixes)) {
            char *new_name = NULL;
            for (int k = 0; k < num_suffixes; k++) {
                char *p = strstr(seq_names[i], suffixes[k]);
                if (p != NULL) {
                    new_name = rename_file(seq_names[i], suffixes[k]);
                    break;
                }
            }
            seq_names[i] = strcpy(seq_names[i], new_name);
            if (new_name != NULL) {
                free(new_name);
            }
        }
    }
    fprintf(fo, "%d\n", num_seqs);
    if (flag == 0) {
        for (int i = 0; i < num_seqs; i++) {
            fprintf(fo, "%s\t", seq_names[i]);
            for (int j = 0; j < num_seqs; j++) {
                if (j <= i) {
                    fprintf(fo, "%.6f\t", distances[i][j]);
                }
            }
            fprintf(fo, "\n");
        }
    } else {
        for (int i = 0; i < num_seqs; i++) {
            fprintf(fo, "%s\t", seq_names[i]);
            for (int j = 0; j < num_seqs; j++) {
                if (j < i) {
                    fprintf(fo, "%.6f\t", distances[i][j]);
                }
            }
            fprintf(fo, "\n");
        }
    }
    for (int i = 0; i < num_seqs; i++) {
        free(distances[i]);
        free(seq_names[i]);
    }
    input_name = NULL;
    output_name = NULL;
    free(seq_names);
    free(distances);
    fclose(fp);
    fclose(fo);
    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("create matrix runtime: %fs\n", cpu_time_used);
    return 1;
}

static PyObject *py_write_dim_shuffle_file(PyObject *self, PyObject *args) {
    int k, s, l;
    char *o;
    if (!PyArg_ParseTuple(args, "iiis", &k, &s, &l, &o)) {
        return NULL;
    }
    dim_shuffle_stat.k = k;
    dim_shuffle_stat.subk = s;
    dim_shuffle_stat.drlevel = l;
    strcpy(shuf_out_file_prefix, o);
    int state = write_dim_shuffle_file(&dim_shuffle_stat, shuf_out_file_prefix);
    return Py_BuildValue("i", state);
}


static PyObject *py_dist_dispatch(PyObject *self, PyObject *args) {
    char *str1;
    char *str2;
    char *str3;
    int flag1;
    int flag2;
    int N;
    if (!PyArg_ParseTuple(args, "sssiii", &str1, &str2, &str3, &flag1, &N, &flag2)) {
        return NULL;
    }
    if (flag1 == 0) {
        dist_opt_val_t dist_opt_val1 =
                {
                        .k = 8,
                        .p = 0,
                        .dr_level = 2,
                        .dr_file = "",
                        .mmry = 0,
                        .fmt = "mfa",
                        .refpath = "",
                        .fpath = "",
                        .outdir = ".",
                        .kmerocrs = 1,
                        .kmerqlty = 0,
                        .keepco = false,
                        .stage2 = false,
                        .num_neigb = 0,
                        .mut_dist_max = 1,
                        .metric = Jcd,
                        .outfields = CI,
                        .correction = false,
                        .abundance = false,
                        .pipecmd = "",
                        .shared_kmerpath="",
                        .keep_shared_kmer=false,
                        .byread=false,
                        .num_remaining_args = 0,
                        .remaining_args = NULL
                };
        struct stat path_stat;
        if (stat(str1, &path_stat) >= 0 && S_ISREG(path_stat.st_mode)) {
            if (strlen(str1) < PATHLEN)
                strcpy(dist_opt_val1.dr_file, str1);
            else
                fprintf(stderr, "-L argument path should not longer than %d", PATHLEN);
        } else {
            if (atoi(str1) >= dist_opt_val1.k - 2 || atoi(str1) < 0)
                fprintf(stderr, "-L: dimension reduction level should never larger than Kmer length - 2,"
                                " which is %d here", dist_opt_val1.k - 2);
            dist_opt_val1.dr_level = atoi(str1);
        }
        strcpy(dist_opt_val1.refpath, str2);
        strcpy(dist_opt_val1.outdir, str3);
        dist_opt_val1.num_remaining_args = 0;
        dist_opt_val1.remaining_args = NULL;
#ifdef _OPENMP
        if(dist_opt_val1.p == 0)
            dist_opt_val1.p = omp_get_num_procs();
#else
        if (dist_opt_val1.p == 0)
            dist_opt_val1.p = 1;
#endif
//#ifdef _WIN32
//        double sys_mm = get_sys_mmry();
//        double rqst_mm = 2.0;
//        if (rqst_mm > sys_mm) {
//            fprintf(stderr, "Memory request is larger than system available %f. Ignoring -m %f", sys_mm, rqst_mm);
//            dist_opt_val1.mmry = sys_mm;
//        } else {
//            dist_opt_val1.mmry = rqst_mm;
//        }
//#endif
        int state1 = dist_dispatch(&dist_opt_val1);
    } else if (flag1 == 1) {
        dist_opt_val_t dist_opt_val2 =
                {
                        .k = 8,
                        .p = 0,
                        .dr_level = 2,
                        .dr_file = "",
                        .mmry = 0,
                        .fmt = "mfa",
                        .refpath = "",
                        .fpath = "",
                        .outdir = ".",
                        .kmerocrs = 1,
                        .kmerqlty = 0,
                        .keepco = false,
                        .stage2 = false,
                        .num_neigb = 0,
                        .mut_dist_max = 1,
                        .metric = Jcd,
                        .outfields = CI,
                        .correction = false,
                        .abundance = false,
                        .pipecmd = "",
                        .shared_kmerpath="",
                        .keep_shared_kmer=false,
                        .byread=false,
                        .num_remaining_args = 0,
                        .remaining_args = NULL
                };
        struct stat path_stat;
        if (stat(str1, &path_stat) >= 0 && S_ISREG(path_stat.st_mode)) {
            if (strlen(str1) < PATHLEN)
                strcpy(dist_opt_val2.dr_file, str1);
            else
                fprintf(stderr, "-L argument path should not longer than %d", PATHLEN);
        }
        strcpy(dist_opt_val2.outdir, str3);
        dist_opt_val2.num_remaining_args = 1;
        dist_opt_val2.remaining_args = &str2;
#ifdef _OPENMP
        if(dist_opt_val2.p == 0)
            dist_opt_val2.p = omp_get_num_procs();
#else
        if (dist_opt_val2.p == 0)
            dist_opt_val2.p = 1;
#endif
//#ifdef _WIN32
//        double sys_mm = get_sys_mmry();
//        double rqst_mm = 2.0;
//        if (rqst_mm > sys_mm) {
//            fprintf(stderr, "Memory request is larger than system available %f. Ignoring -m %f", sys_mm, rqst_mm);
//            dist_opt_val2.mmry = sys_mm;
//        } else {
//            dist_opt_val2.mmry = rqst_mm;
//        }
//#endif
        int state2 = dist_dispatch(&dist_opt_val2);
    } else if (flag1 == 2) {
        dist_opt_val_t dist_opt_val3 =
                {
                        .k = 8,
                        .p = 0,
                        .dr_level = 2,
                        .dr_file = "",
                        .mmry = 0,
                        .fmt = "mfa",
                        .refpath = "",
                        .fpath = "",
                        .outdir = ".",
                        .kmerocrs = 1,
                        .kmerqlty = 0,
                        .keepco = false,
                        .stage2 = false,
                        .num_neigb = 0,
                        .mut_dist_max = 1,
                        .metric = Jcd,
                        .outfields = CI,
                        .correction = false,
                        .abundance = false,
                        .pipecmd = "",
                        .shared_kmerpath="",
                        .keep_shared_kmer=false,
                        .byread=false,
                        .num_remaining_args = 0,
                        .remaining_args = NULL
                };
        strcpy(dist_opt_val3.refpath, str1);
        char *output_name = "distout";
        strcpy(dist_opt_val3.outdir, output_name);
        dist_opt_val3.num_remaining_args = 1;
        dist_opt_val3.remaining_args = &str3;
        dist_opt_val3.num_neigb = N;
#ifdef _OPENMP
        if(dist_opt_val3.p == 0)
            dist_opt_val3.p = omp_get_num_procs();
#else
        if (dist_opt_val3.p == 0)
            dist_opt_val3.p = 1;
#endif
//#ifdef _WIN32
//        double sys_mm = get_sys_mmry();
//        double rqst_mm = 6.5;
//        printf("sys_mm: %f GB\n", sys_mm);
//        printf("rqst_mm: %f GB\n", rqst_mm);
//        if (rqst_mm > sys_mm) {
//            fprintf(stderr, "Memory request is larger than system available %f. Ignoring -m %f", sys_mm, rqst_mm);
//            dist_opt_val3.mmry = sys_mm;
//        } else {
//            dist_opt_val3.mmry = rqst_mm;
//        }
//#endif
        int state3 = dist_dispatch(&dist_opt_val3);
        if (state3 == 0) {
            if (flag2 != 3) {
                char *input_name = "distout/distance.out";
                char cwd[1024];
                if (getcwd(cwd, sizeof(cwd)) != NULL) {
                    char final_path[1024];
                    strcpy(final_path, cwd);
                    strcat(final_path, "/");
                    strcat(final_path, input_name);
                    create_matrix(final_path, str2, flag2);
                }
            }
        }
    } else {
        dist_opt_val_t dist_opt_val4 =
                {
                        .k = 8,
                        .p = 0,
                        .dr_level = 2,
                        .dr_file = "",
                        .mmry = 0,
                        .fmt = "mfa",
                        .refpath = "",
                        .fpath = "",
                        .outdir = ".",
                        .kmerocrs = 1,
                        .kmerqlty = 0,
                        .keepco = false,
                        .stage2 = false,
                        .num_neigb = 0,
                        .mut_dist_max = 1,
                        .metric = Jcd,
                        .outfields = CI,
                        .correction = false,
                        .abundance = false,
                        .pipecmd = "",
                        .shared_kmerpath="",
                        .keep_shared_kmer=false,
                        .byread=false,
                        .num_remaining_args = 0,
                        .remaining_args = malloc(2 * sizeof(char *))
                };
        strcpy(dist_opt_val4.outdir, str1);
        dist_opt_val4.num_remaining_args = 2;
        dist_opt_val4.remaining_args[0] = str2;
        dist_opt_val4.remaining_args[1] = str3;
#ifdef _OPENMP
        if(dist_opt_val4.p == 0)
            dist_opt_val4.p = omp_get_num_procs();
#else
        if (dist_opt_val4.p == 0)
            dist_opt_val4.p = 1;
#endif
//#ifdef _WIN32
//        double sys_mm = get_sys_mmry();
//        double rqst_mm = 2.0;
//        if (rqst_mm > sys_mm) {
//            fprintf(stderr, "Memory request is larger than system available %f. Ignoring -m %f", sys_mm, rqst_mm);
//            dist_opt_val4.mmry = sys_mm;
//        } else {
//            dist_opt_val4.mmry = rqst_mm;
//        }
//#endif
        int state4 = dist_dispatch(&dist_opt_val4);
    }
    str1 = NULL;
    str2 = NULL;
    str3 = NULL;
    return Py_BuildValue("i", 1);
}


static PyObject *py_sketch_union(PyObject *self, PyObject *args) {
    char *i;
    char *o;
    if (!PyArg_ParseTuple(args, "ss", &i, &o)) {
        return NULL;
    }
    set_opt_t set_opt1 = {
            .operation = -1,
            .p = 1,
            .P = 0,
            .num_remaining_args = 0,
            .remaining_args = NULL,
            .insketchpath = "",
            .pansketchpath="",
            .subsetf[0] = '\0',
            .outdir = ""
    };
    strcpy(set_opt1.insketchpath, i);
    strcpy(set_opt1.outdir, o);
    int state = sketch_union(set_opt1);
    i = NULL;
    o = NULL;
    return Py_BuildValue("i", state);
}


static PyObject *py_sketch_sub(PyObject *self, PyObject *args) {
    char *i;
    char *o;
    char *remain;
    if (!PyArg_ParseTuple(args, "sss", &i, &o, &remain)) {
        return NULL;
    }
    set_opt_t set_opt2 = {
            .operation = -1,
            .p = 1,
            .P = 0,
            .num_remaining_args = 0,
            .remaining_args = NULL,
            .insketchpath = "",
            .pansketchpath="",
            .subsetf[0] = '\0',
            .outdir = ""
    };
    strcpy(set_opt2.pansketchpath, i);
    strcpy(set_opt2.outdir, o);
    strcpy(set_opt2.insketchpath, remain);
    set_opt2.num_remaining_args = 1;
    set_opt2.remaining_args = &remain;
    int state = sketch_operate(set_opt2);
    set_opt2.num_remaining_args = 0;
    set_opt2.remaining_args = NULL;
    i = NULL;
    o = NULL;
    remain = NULL;
    return Py_BuildValue("i", state);
}


static PyObject *py_sketch_print(PyObject *self, PyObject *args) {
    char *o;
    char *remain;
    if (!PyArg_ParseTuple(args, "ss", &remain, &o)) {
        return NULL;
    }
    set_opt_t set_opt3 = {
            .operation = -1,
            .p = 1,
            .P = 0,
            .num_remaining_args = 0,
            .remaining_args = NULL,
            .insketchpath = "",
            .pansketchpath="",
            .subsetf[0] = '\0',
            .outdir = ""
    };
    set_opt3.P = 1;
    strcpy(set_opt3.insketchpath, remain);
    set_opt3.num_remaining_args = 1;
    set_opt3.remaining_args = &remain;
    print_gnames(set_opt3, o);
    set_opt3.num_remaining_args = 0;
    set_opt3.remaining_args = NULL;
    o = NULL;
    remain = NULL;
    return Py_BuildValue("i", 1);
}


static PyObject *py_sketch_group(PyObject *self, PyObject *args) {
    char *i;
    char *o;
    char *remain;
    if (!PyArg_ParseTuple(args, "sss", &i, &remain, &o)) {
        return NULL;
    }
    set_opt_t set_opt4 = {
            .operation = -1,
            .p = 1,
            .P = 0,
            .num_remaining_args = 0,
            .remaining_args = NULL,
            .insketchpath = "",
            .pansketchpath="",
            .subsetf[0] = '\0',
            .outdir = ""
    };
    strcpy(set_opt4.subsetf, i);
    strcpy(set_opt4.outdir, o);
    strcpy(set_opt4.insketchpath, remain);
    set_opt4.num_remaining_args = 1;
    set_opt4.remaining_args = &remain;
    int state = grouping_genomes(set_opt4, i);
    set_opt4.num_remaining_args = 0;
    set_opt4.remaining_args = NULL;
    i = NULL;
    o = NULL;
    remain = NULL;
    return Py_BuildValue("i", state);
}


static PyMethodDef KssdMethods[] = {
        {"write_dim_shuffle_file", py_write_dim_shuffle_file, METH_VARARGS, "shuffle"},
        {"dist_dispatch",          py_dist_dispatch,          METH_VARARGS, "sketch and dist"},
        {"sketch_union",           py_sketch_union,           METH_VARARGS, "sketch union"},
        {"sketch_operate",         py_sketch_sub,             METH_VARARGS, "sketch sub"},
        {"print_gnames",           py_sketch_print,           METH_VARARGS, "sketch print"},
        {"grouping_genomes",       py_sketch_group,           METH_VARARGS, "sketch group"},
        {NULL, NULL,                                          0, NULL}
};

static struct PyModuleDef kssdmodule = {
        PyModuleDef_HEAD_INIT,
        "kssd",
        "A kssd module",
        -1,
        KssdMethods
};

PyMODINIT_FUNC PyInit_kssd(void) {
    return PyModule_Create(&kssdmodule);
}