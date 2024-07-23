#include <Python.h>
#include "njheaders/util.h"
#include "njheaders/align.h"
#include "njheaders/cluster.h"
#include "njheaders/buildtree.h"
#include "njheaders/distancemat.h"
#include "njheaders/tree.h"


void build(FILE *input, const char *output_name) {
    struct DistanceMatrix *mat;
    struct Alignment *aln;
    struct Tree *myTree;
    struct ClusterGroup *group;
    //step1
    mat = read_phylip_DistanceMatrix(input, &aln);
    fclose(input);
    group = alignment_to_ClusterGroup(aln, FALSE);
    group->matrix = mat;
    //step2
    myTree = neighbour_joining_buildtree(group, FALSE);
    FILE *handle = fopen(output_name, "w");
    if (handle == NULL) {
        printf("Failed to open file\n");
        return;
    } else {
        write_newhampshire_Tree(handle, myTree, FALSE);
        fclose(handle);
    }
    aln = free_Alignment(aln);
    group = free_ClusterGroup(group);
    myTree = free_Tree(myTree);
}

static PyObject *py_build(PyObject *self, PyObject *args) {
    char *input_name;
    char *output_name;
    FILE *matrixfile;
    if (!PyArg_ParseTuple(args, "ss", &input_name, &output_name)) {
        return NULL;
    }
    matrixfile = fopen(input_name, "r");
    if (matrixfile == NULL)
        fatal_util("Could not open file %s for reading", input_name);
    build(matrixfile, output_name);
    return Py_BuildValue("i", 1);
}

static PyMethodDef NJMethods[] = {
        {"build", py_build, METH_VARARGS, "build"},
        {NULL, NULL,                0, NULL}
};

static struct PyModuleDef njmodule = {
        PyModuleDef_HEAD_INIT,
        "nj",
        "A nj module",
        -1,
        NJMethods
};

PyMODINIT_FUNC PyInit_nj(void) {
    return PyModule_Create(&njmodule);
}