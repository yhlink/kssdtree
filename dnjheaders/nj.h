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

#ifndef NJ
typedef struct njThread NJthread;
struct njThread {
	pthread_t id;
	Matrix *D;
	Vector *sD;
	Vector *Q;
	int *N;
	int *P;
	double min;
	int mi;
	int mj;
	long unsigned pos;
	unsigned num;
	struct njThread *next;
};
#define NJ 1
#endif

extern void (*limbLengthPtr)(double*, double*, int, int, Vector *, int*, double);
extern double (*initQchunkPtr)(Matrix *, double*, int*, double, int*, int*, int, int);
extern double (*minDchunkPtr)(Matrix *, int*, double, int*, int*, int, int);
extern void (*updateDptr)(Matrix *, Vector *, int *, int, int, double, double);
extern long unsigned (*minDist)(Matrix *, Vector *, int *);
extern void * (*minDist_thread)(void *);

void limbLength(double *Li, double *Lj, int i, int j, Vector *sD, int *N, double D_ij);
void limbLengthNeg(double *Li, double *Lj, int i, int j, Vector *sD, int *N, double D_ij);
int * initSummaD(Vector *sD, Matrix *D, int *N);
long unsigned initQ(Matrix *D, Vector *sD, int *N);
double initQchunk(Matrix *D, double *sD, int *N, double min, int *mi, int *mj, int i, int j);
long unsigned initQ_MN(Matrix *D, Vector *sD, int *N);
double initQ_MNchunk(Matrix *D, double *sD, int *N, double max, int *mi, int *mj, int i, int j);
void * initQ_thread(void *arg);
long unsigned minD(Matrix *D, Vector *sD, int *N);
double minDchunk(Matrix *D, int *N, double min, int *mi, int *mj, int i, int j);
long unsigned maxD(Matrix *D, Vector *sD, int *N);
double maxDchunk(Matrix *D, int *N, double max, int *mi, int *mj, int i, int j);
void * minD_thread(void *arg);
long unsigned minPair(Matrix *Q);
void updateD(Matrix *D, Vector *sD, int *N, int i, int j, double Li, double Lj);
void updateD_UPGMA(Matrix *D, Vector *sD, int *N, int i, int j, double Li, double Lj);
void updateD_FF(Matrix *D, Vector *sD, int *N, int i, int j, double Li, double Lj);
void updateD_CF(Matrix *D, Vector *sD, int *N, int i, int j, double Li, double Lj);
int * nj(Matrix *D, Vector *sD, int *N, Qseqs **names);
int * nj_thread(Matrix *D, Vector *sD, int *N, Qseqs **names, int thread_num);
