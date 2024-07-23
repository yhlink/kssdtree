/* Philip T.L.C. Clausen May 2021 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "matrix.h"
#include "qseqs.h"
#include "vector.h"

extern long unsigned (*Qpair)(Matrix *, Vector *, Vector *, int *, int *, int);
extern long unsigned (*Qrow)(Matrix *, Vector *, Vector *, int *, int *, double *, int);
extern int (*Qbool)(double, double, long unsigned, long unsigned);
extern int (*nextQrow)(int, double, double*);
extern int (*qPos)(double *, int, int);

long unsigned minQpair(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i);
long unsigned maxQpair(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i);
long unsigned UPGMApair(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i);
int nextQminRow(int i, double min, double *Q);
int nextQmaxRow(int i, double max, double *Q);
long unsigned minQrow(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, double *min, int i);
long unsigned maxQrow(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, double *max, int i);
long unsigned UPGMArow(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, double *min, int i);
int minQbool(double min, double Min, long unsigned pos, long unsigned Pos);
int maxQbool(double max, double Max, long unsigned pos, long unsigned Pos);
void * minQ_thread(void *arg);
int updateDNJ(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj);
int updateDMN(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj);
int DNJ_popArrange(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int pos);
int minPos(double *Q, int i, int j);
int maxPos(double *Q, int i, int j);
int * dnj(Matrix *D, Vector *sD, Vector *Q, int *N, Qseqs **names);
int * dnj_thread(Matrix *D, Vector *sD, Vector *Q, int *N, Qseqs **names, int thread_num);
