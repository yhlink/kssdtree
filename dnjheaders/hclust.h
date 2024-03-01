/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

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

#include <pthread.h>
#include "matrix.h"
#include "qseqs.h"
#include "vector.h"

extern int * (*initDsDQN)(Matrix *, Vector *, Vector *, int *);
extern long unsigned (*pairQ)(Vector *, int *);
extern int (*updateDsDQNPtr)(Matrix *, Vector *, Vector *, int *, int *, int, int, double, double);
extern int (*popArrangePtr)(Matrix *, Vector *, Vector *, int*, int*, int);

int * initAlloc(Matrix *D, Vector *sD, Vector *Q, int *N);
int * initHNJ(Matrix *D, Vector *sD, Vector *Q, int *N);
int * initHMN(Matrix *D, Vector *sD, Vector *Q, int *N);
int * initDmin(Matrix *D, Vector *sD, Vector *Q, int *N);
int * initDmax(Matrix *D, Vector *sD, Vector *Q, int *N);
long unsigned minQ(Vector *Q, int *P);
long unsigned maxQ(Vector *Q, int *P);
void updatePrevQ(Matrix *D, Vector *sD, Vector *Q, int *N, int *P);
int updateHNJ(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj);
int updateHMN(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj);
int updateUPGMA(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj);
int updateFF(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj);
int updateCF(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj);
int HNJ_popArrange(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int pos);
int HMN_popArrange(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int pos);
int UPGMA_popArrange(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int pos);
int * hclust(Matrix *D, Vector *sD, Vector *Q, int *N, Qseqs **names);
