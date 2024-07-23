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
#include "kssdheaders/command_dist_wrapper.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#endif

const char *mk_dist_rslt_dir(const char *parentdirpath, const char *outdirpath) {
    struct stat dstat;
    const char *outfullpath = malloc(PATHLEN * sizeof(char));
    sprintf((char *) outfullpath, "%s/%s", parentdirpath, outdirpath);
    if (stat(parentdirpath, &dstat) == 0 && S_ISDIR(dstat.st_mode)) {
        if (stat(outfullpath, &dstat) == 0) {
            errno = EEXIST;
            fprintf(stderr, "%s", outfullpath);
        } else {
#ifdef _WIN32
            mkdir(outfullpath);
#else
            mkdir(outfullpath, 0777);
#endif
        }
    } else {
#ifdef _WIN32
        mkdir(parentdirpath);
        mkdir(outfullpath);
#else
        mkdir(parentdirpath, 0777);
        mkdir(outfullpath, 0777);
#endif
    }
    return outfullpath;
};
