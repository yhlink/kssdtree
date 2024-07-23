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
#include "kssdheaders/command_shuffle.h"
#include "kssdheaders/global_basic.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>

#ifdef _WIN32
void init_genrand(uint32_t seed) {
    msstate[0] = seed;
    for (unsigned int i = 1; i < 624; ++i) {
        msstate[i] = 0x6c078965 * (msstate[i - 1] ^ (msstate[i - 1] >> 30)) + i;
    }
}

uint32_t genrand() {
    if (msindex == 0) {
        for (unsigned int i = 0; i < 624; ++i) {
            uint32_t y = (msstate[i] & 0x80000000) + (msstate[(i + 1) % 624] & 0x7fffffff);
            msstate[i] = msstate[(i + 397) % 624] ^ (y >> 1);
            if (y % 2 != 0) {
                msstate[i] ^= 0x9908b0df;
            }
        }
    }
    uint32_t y = msstate[msindex];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= (y >> 18);
    msindex = (msindex + 1) % 624;
    return y;
}
#endif



int *shuffle(int arr[], int len_arr) {
    if (len_arr > 2147483647)
        fprintf(stderr, "shuffling array length %d longer than RAND_MAX: %d", len_arr, RAND_MAX);
#ifdef _WIN32
    init_genrand((uint32_t) time(0));
    int j, temp;
    for (int i = len_arr - 1; i > 0; i--) {
        j = genrand() % (i + 1);
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
#else
//    printf("Linux/MacOS\n");
    srand(time(NULL));
    int j, temp;
    for (int i = len_arr - 1; i > 0; i--) {
        j = rand() % (i + 1);
        temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }
#endif
    return arr;
};

int *shuffleN(int n, int base) {
    int *arr;
    arr = malloc(n * sizeof(int));
    if (arr == NULL) fprintf(stderr, "shuffleN");
    for (int i = 0; i < n; i++)
        arr[i] = i + base;
    return shuffle(arr, n);
};

int add_len_drlevel2subk(void) {
    int min_smp_len = 0, min_subctx_dim_smp_sz;
    min_subctx_dim_smp_sz = MIN_SUBCTX_DIM_SMP_SZ;
    while (min_subctx_dim_smp_sz >>= 1) { min_smp_len++; };
    return ceil((float) min_smp_len / 4);
};

int write_dim_shuffle_file(dim_shuffle_stat_t *dim_shuffle_stat, char *outfile_prefix) {
    if (dim_shuffle_stat->k < dim_shuffle_stat->subk)
        fprintf(stderr, "write_dim_shuffle_file(): half-context len: %d"
                   "should larger than half-subcontext len (or dimension reduce level + 2) %d",
            dim_shuffle_stat->k, dim_shuffle_stat->subk);
    if (dim_shuffle_stat->subk >= 8)
        fprintf(stderr, "write_dim_shuffle_file(): subk shoud smaller than 8");
    int dim_after_reduction = 1 << 4 * (dim_shuffle_stat->subk - dim_shuffle_stat->drlevel);
    if (dim_after_reduction < MIN_SUBCTX_DIM_SMP_SZ)
        printf("dimension after reduction %d is smaller than the suggested minimal"
              " dimension sample size %d, which might cause loss of robutness, -s %d is suggested\n", dim_after_reduction,
              MIN_SUBCTX_DIM_SMP_SZ, dim_shuffle_stat->drlevel + 3);
    if (strlen(outfile_prefix) + strlen(".shuf") > PATHLEN)
        fprintf(stderr, "output path name %s should less than %d characters", outfile_prefix, PATHLEN);
    char outfile[PATHLEN];
    sprintf(outfile, "%s.shuf", outfile_prefix);
    FILE *shuf_out;
    if ((shuf_out = fopen(outfile, "wb")) == NULL)
        fprintf(stderr, "write_dim_shuffle_file(): open file %s failed", outfile);
#ifdef _WIN32
    init_genrand((uint32_t) time(0));
    dim_shuffle_stat->id = genrand();
#else
    srand(time(NULL));
    dim_shuffle_stat->id = rand();
#endif
    int *shuffled_dim = shuffleN(1 << 4 * dim_shuffle_stat->subk, 0);
    shuffled_dim = shuffle(shuffled_dim, 1 << 4 * dim_shuffle_stat->subk);
    int ret = fwrite(dim_shuffle_stat, sizeof(dim_shuffle_stat_t), 1, shuf_out)
              + fwrite(shuffled_dim, sizeof(int), 1 << 4 * dim_shuffle_stat->subk, shuf_out);
    fclose(shuf_out);
    free(shuffled_dim);
    return ret;
};

dim_shuffle_t *read_dim_shuffle_file(char *dim_shuffle_file) {
    int basename_len = strlen(dim_shuffle_file) - strlen(".shuf");
    if (strcmp((dim_shuffle_file + basename_len), ".shuf") != 0)
        fprintf(stderr, "read_dim_shuffle_file(): input file %s is not .shuf file", dim_shuffle_file);
    FILE *shuf_in;
    if ((shuf_in = fopen(dim_shuffle_file, "rb")) == NULL)
        fprintf(stderr, "read_dim_shuffle_file(): open file %s failed", dim_shuffle_file);
    dim_shuffle_t *dim_shuffle = malloc(sizeof(dim_shuffle_t));
    fread(&(dim_shuffle->dim_shuffle_stat), sizeof(dim_shuffle_stat_t), 1, shuf_in);
    int shuf_arr_len = 1 << 4 * dim_shuffle->dim_shuffle_stat.subk;
    dim_shuffle->shuffled_dim = malloc(sizeof(int) * shuf_arr_len);
    fread(dim_shuffle->shuffled_dim, sizeof(int) * shuf_arr_len, 1, shuf_in);
    fclose(shuf_in);
    return dim_shuffle;
}
