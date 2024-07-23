#ifndef TREE_COMMAND_SET_H
#define TREE_COMMAND_SET_H
#include "global_basic.h"
typedef struct subset {
    int taxid;
    char *taxname;
    int *gids;
} subset_t;
typedef struct compan {
    int taxn;
    int gn;
    subset_t *tax;
} compan_t;
typedef struct set_opt {
    int operation;
    int p;
    int P;
    int num_remaining_args;
    char **remaining_args;
    char insketchpath[PATHLEN];
    char pansketchpath[PATHLEN];
    char subsetf[PATHLEN];
    char outdir[PATHLEN];
} set_opt_t;
int sketch_union(set_opt_t set_opt);

int sketch_operate(set_opt_t set_opt);

int uniq_sketch_union(set_opt_t set_opt);

int combin_pans(set_opt_t set_opt);

void print_gnames(set_opt_t set_opt, char *outputname);

compan_t *organize_taxf(char* taxfile);
int combin_subset_pans(set_opt_t set_opt,char* taxfile);
int grouping_genomes(set_opt_t set_opt,char* taxfile);
#endif //TREE_COMMAND_SET_H
