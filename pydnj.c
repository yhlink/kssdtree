#include <Python.h>
#include <stdio.h>
#include <time.h>
#include "dnjheaders/bytescale.h"
#include "dnjheaders/dnj.h"
#include "dnjheaders/filebuff.h"
#include "dnjheaders/hclust.h"
#include "dnjheaders/matrix.h"
#include "dnjheaders/nj.h"
#include "dnjheaders/nwck.h"
#include "dnjheaders/pherror.h"
#include "dnjheaders/phy.h"
#include "dnjheaders/qseqs.h"
#include "dnjheaders/tmp.h"
#include "dnjheaders/vector.h"

//#ifdef _OPENMP
//#include <omp.h>
//#else
//#define omp_get_thread_num() 0
//#endif

void build(char *inputfilename, char *outputfilename, int flag, char sep, char quotes, char m, int thread_num) {
//    printf("thread_num: %d\n", thread_num);
    int i, *N;
    FILE *outfile;
    FileBuff *infile;
    Matrix *D;
    Qseqs **names, *header;
    Vector *sD, *Q;
    time_t t0, t1;

    /* init */
    outfile = (*outputfilename == '-' && outputfilename[1] == 0) ? stdout : sfopen(outputfilename, "wb");
    infile = setFileBuff(1048576);
    header = setQseqs(64);
    i = 32;
    D = ltdMatrix_init(i);
    sD = vector_init(i);
    if (m == 'e') {
        Q = 0;
        N = smalloc(i * sizeof(int));
    } else {
        Q = vector_init(i);
        N = smalloc(2 * i * sizeof(int));
    }
    names = smalloc(i * sizeof(Qseqs *));
    names += i;
    ++i;
    while (--i) {
        *--names = setQseqs(4);
    }
    /* set ptr according ot flag */
    if (flag) {
        if (flag & 1) {
            formLastNodePtr = &formLastBiNode;
        }
        if (flag & 2) {
            limbLengthPtr = &limbLengthNeg;
        }
    }
    /* set */
    openAndDetermine(infile, inputfilename);

    /* generate trees */
    t0 = clock();
    while ((names = loadPhy(D, names, header, infile, sep, quotes)) && D->n) {
        t1 = clock();
        fprintf(stderr, "# Total time used loading matrix: %.2f s.\n", difftime(t1, t0) / 1000000);
        t0 = t1;
        if (2 < D->n) {
            /* make tree */
            if (m == 'd') {
                N = dnj_thread(D, sD, Q, N, names, thread_num);
            } else if (m == 'e') {
                N = nj_thread(D, sD, N, names, thread_num);
            } else { /* m == 'h' */
                N = hclust(D, sD, Q, N, names);
            }
        } else if (D->n == 2) {
            /* form tree */
            formLastBiNode(*names, names[1],
                           (D->mat ? **(D->mat) : D->fmat ? **(D->fmat) : D->smat ? uctod(**(D->smat)) : uctod(
                                   **(D->bmat))));
        }

        /* output tree */
        if (header->len) {
            fprintf(outfile, ">%s%s;\n", header->seq, (*names)->seq);
        } else {
            fprintf(outfile, "%s;\n", (*names)->seq);
        }

        t1 = clock();
        fprintf(stderr, "# Total time used Constructing tree: %.2f s.\n", difftime(t1, t0) / 1000000);
        t0 = t1;
    }
    /* clean */
    fclose(outfile);
    closeFileBuff(infile);
    free(N);
    Matrix_destroy(D);
    vector_destroy(sD);
    free(names);
    destroyFileBuff(infile);
}

static PyObject *py_build(PyObject *self, PyObject *args) {
    char *inputfilename;
    char *outputfilename;
    char *method;
    FILE *file;
    char m, sep, quotes;
    int flag, p;
    sep = '\t';
    quotes = '\0';
    flag = 0;
    p = 1;
    int size, precision;
    size = sizeof(double);
    precision = 9;
    if (!PyArg_ParseTuple(args, "sss", &inputfilename, &outputfilename, &method)) {
        return NULL;
    }
    file = fopen(inputfilename, "r");
    if (file == NULL)
        fprintf(stderr, "Could not open file %s for reading", inputfilename);
    /* set print precision */
    setPrecisionNwck(precision);
    /* set precision */
    ltdMatrixInit(-size);
    ltdMatrixMinit(-size);
    if (strcmp(method, "hnj") == 0) {
        /* heuristic neighbor-joining */
//        printf("s0-hnj:%c\n", m);
        m = 'h';
        initDsDQN = &initHNJ;
        pairQ = &minQ;
        updateDsDQNPtr = &updateHNJ;
        popArrangePtr = &HNJ_popArrange;
    } else if (strcmp(method, "dnj") == 0) {
        /* dynamic neighbor-joining */
//        printf("s0-dnj:%c\n", m);
        m = 'd';
        Qpair = &minQpair;
        minDist_thread = &minQ_thread;
        Qrow = &minQrow;
        nextQrow = &nextQminRow;
        Qbool = &minQbool;
        initDsDQN = &initHNJ;
        pairQ = &minQ;
        updateDsDQNPtr = &updateDNJ;
        popArrangePtr = &DNJ_popArrange;
        qPos = &minPos;
    } else {
        /* neighbor-joining */
//        printf("s0-nj:%c\n", m);
        m = 'e';
        updateDptr = &updateD;
        minDist = &initQ;
        minDist_thread = &initQ_thread;
        initQchunkPtr = &initQchunk;
    }
#ifdef _OPENMP
    p = omp_get_num_procs();
#endif
    printf("p=%d\n", p);
    build(inputfilename, outputfilename, flag, sep, quotes, m, p);
    free(file);
    return Py_BuildValue("i", 1);
}

static PyMethodDef DNJMethods[] = {
        {"build", py_build, METH_VARARGS, "build"},
        {NULL, NULL,        0, NULL}
};

static struct PyModuleDef dnjmodule = {
        PyModuleDef_HEAD_INIT,
        "dnj",
        "A dnj module",
        -1,
        DNJMethods
};

PyMODINIT_FUNC PyInit_dnj(void) {
    return PyModule_Create(&dnjmodule);
}