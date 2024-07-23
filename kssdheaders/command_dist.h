//	Copyright 2019 Huiguang Yi. All Rights Reservered.
//
//	Licensed under the Apache License, Version 2.0 (the "License");
//	you may not use this file except in compliance with the License.
//	You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//	Unless required by applicable law or agreed to in writing, software
//	distributed under the License is distributed on an "AS IS" BASIS,
//	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//	See the License for the specific language governing permissions and
//	limitations under the License.
#ifndef COMMAND_DIST
#define COMMAND_DIST

#include "global_basic.h"
#include "command_dist_wrapper.h"
#include "command_shuffle.h"
#include <string.h>

#define DISM_MEM_PCT 0.25
#define CO_MEM_PCT 0.125
extern const char co_dstat[];
typedef struct mem_dispatch {
    bool has_ref;
    bool has_arg;
    bool dm_inmem;
    bool keep_coindisk;
    bool keep_mcoinmem;
    bool alphabet_group;
    bool multiple_composition;
    enum {
        OBJ_COMPATIBLE = 4096,
        OBJ_TRADITION = 16384,
        OBJ_GROUP = 65536,
    } enum_bin_sz;
    enum {
        SLIM,
        STD,
        ALL,
    } gbin;
    enum {
        MAKE_REF,
        REF_QRY,
        DIST,
    } usage_mode;
    bool mco_info;
} mem_dispatch_t;
typedef struct mem_usage_stat {
    llong shuffled_subctx_arr_sz;
    llong input_file_name_sz;
    llong others;
} mem_usage_stat_t;
typedef long long int llint;
typedef struct mco_dirstat {
    unsigned int shuf_id;
    int kmerlen;
    int dim_rd_len;
    int comp_num;
    int infile_num;
} mco_dstat_t;
typedef struct output {
    int X_size;
    int Y_size;
    int XnY_size;
    double metric;
    double dist;
    double CI95_mtrc_1;
    double CI95_mtrc_2;
    double CI95_dist_1;
    double CI95_dist_2;
    double pv;
} output_t;
typedef struct print_ctrl {
    MTRIC metric;
    PFIELD pfield;
    bool correction;
    double dthreshold;
    unsigned int Y_size;
    llong cmprsn_num;
    char *qname;
    int qry_len;
} print_ctrl_t;
extern dim_shuffle_t *dim_shuffle;
extern unsigned int hashsize;
extern int component_num;

int dist_dispatch(struct dist_opt_val *opt_val);

infile_tab_t *dist_organize_infiles(dist_opt_val_t *opt_val);

infile_tab_t *dist_organize_refpath(dist_opt_val_t *opt_val);

dim_shuffle_t *get_dim_shuffle(dist_opt_val_t *opt_val_in);

int get_hashsz(dim_shuffle_t *dim_shuffle_in);

const char *run_stageI(dist_opt_val_t *opt_val, infile_tab_t *seqfile_stat,
                       int *shuffled_seqfname_ind, const char *co_dir, int p_fit_mem);

void run_stageII(const char *co_dstat_fpath, const char *dist_mco_dir, int p_fit_mem);

void mco_cbdco_nobin_dist(dist_opt_val_t *opt_val_in);

void mco_co_dist(char *refmco_dname, char *qryco_dname, const char *distout_dir, int p_fit_mem);

void mco_cbd_co_dist(dist_opt_val_t *opt_val_in);

void mco_cbd_koc_compatible_dist(dist_opt_val_t *opt_val_in);

void dist_print(const char *distf, FILE *dist_fp);

void fname_dist_print(int ref_bin_code, int qry_fcode, const char *distout_dir, unsigned int *ref_ctx_ct_list,
                      unsigned int *qry_ctx_ct_list, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN],
                      FILE *dout_fp);

void
dist_print_nobin(const char *distout_dir, unsigned int ref_num, unsigned int qry_num, unsigned int *ref_ctx_ct_list,
                 unsigned int *qry_ctx_ct_list, int num_cof_batch, char (*refname)[PATHLEN], char (*qryfname)[PATHLEN],
                 dist_opt_val_t *opt_val);

void
koc_dist_print_nobin(const char *distout_dir, unsigned int ref_num, unsigned int qry_num, unsigned int *ref_ctx_ct_list,
                     unsigned int *qry_ctx_ct_list, int num_cof_batch, char (*refname)[PATHLEN],
                     char (*qryfname)[PATHLEN]);

const char *combine_queries(dist_opt_val_t *opt_val);

#endif
