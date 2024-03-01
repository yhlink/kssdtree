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

#define _XOPEN_SOURCE 600
#include <float.h>
#include <limits.h>
#include <pthread.h>
#include <stdlib.h>
#include "dnjheaders/bytescale.h"
#include "dnjheaders/matrix.h"
#include "dnjheaders/nj.h"
#include "dnjheaders/nwck.h"
#include "dnjheaders/qseqs.h"
#include "dnjheaders/pherror.h"
#include "dnjheaders/str.h"
#include "dnjheaders/threader.h"
#include "dnjheaders/vector.h"

void (*limbLengthPtr)(double*, double*, int, int, Vector *, int*, double) = &limbLength;
double (*initQchunkPtr)(Matrix *, double*, int*, double, int*, int*, int, int) = &initQchunk;
double (*minDchunkPtr)(Matrix *, int*, double, int*, int*, int, int) = &minDchunk;
void (*updateDptr)(Matrix *, Vector *, int*, int, int, double, double) = &updateD;
long unsigned (*minDist)(Matrix *, Vector *, int*) = &initQ;
void * (*minDist_thread)(void *) = &initQ_thread;

void limbLength(double *Li, double *Lj, int i, int j, Vector *sD, int *N, double D_ij) {
	
	int Ni, Nj;
	double delta_ij;
	
	Ni = N[i] - 2;
	Nj = N[j] - 2;
	if(0 < Ni && 0 < Nj) {
		/* semi-standard */
		delta_ij = ((sD->vec[i] - D_ij) / Ni) - ((sD->vec[j] - D_ij) / Nj);
		*Li = (D_ij + delta_ij) / 2;
		*Lj = (D_ij - delta_ij) / 2;
		
		/* org:
		delta_ij = (sD->vec[i] - sD->vec[j]) / N = (sD->vec[i] - D_ij) / N) - ((sD->vec[j] - D_ij) / N;
		*/
		
		/* avoid negative branch lengths */
		if(*Li < 0) {
			*Lj = D_ij;
			*Li = 0;
		} else if(*Lj < 0) {
			*Li = D_ij;
			*Lj = 0;
		}
	} else if(0 < Ni) {
		/* j only shares similarity with i */
		*Li = 0;
		*Lj = D_ij;
	} else if(0 < Nj) {
		/* i only shares similarity with j */
		*Li = D_ij;
		*Lj = 0;
	} else {
		/* i and j only share similarity with eachother */
		*Li = (*Lj = D_ij / 2);
	}
}

void limbLengthNeg(double *Li, double *Lj, int i, int j, Vector *sD, int *N, double D_ij) {
	
	int Ni, Nj;
	double delta_ij;
	
	Ni = N[i] - 2;
	Nj = N[j] - 2;
	if(0 < Ni && 0 < Nj) {
		/* semi-standard */
		delta_ij = ((sD->vec[i] - D_ij) / Ni) - ((sD->vec[j] - D_ij) / Nj);
		*Li = (D_ij + delta_ij) / 2;
		*Lj = (D_ij - delta_ij) / 2;
		
		/* org:
		delta_ij = (sD->vec[i] - sD->vec[j]) / N = (sD->vec[i] - D_ij) / N) - ((sD->vec[j] - D_ij) / N;
		*/
	} else if(0 < Ni) {
		/* j only shares similarity with i */
		*Li = 0;
		*Lj = D_ij;
	} else if(0 < Nj) {
		/* i only shares similarity with j */
		*Li = D_ij;
		*Lj = 0;
	} else {
		/* i and j only share similarity with eachother */
		*Li = (*Lj = D_ij / 2);
	}
}

int * initSummaD(Vector *sD, Matrix *D, int *N) {
	
	int i, j, *Nptr, *Nvec;
	double dist, *sDptr, *sDvec, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/*
	sD(i) = summa(d(i,k));
	*/
	
	/* alloc */
	if(sD->size < (sD->n = D->n)) {
		vector_realloc(sD, sD->n);
		if(!(N = realloc(N, sD->n * sizeof(int)))) {
				ERROR();
		}
	}
	
	/* init to zero */
	i = sD->n + 1;
	sDptr = sD->vec - 1;
	Nptr = N - 1;
	while(--i) {
		*++sDptr = 0;
		*++Nptr = 1;
	}
	
	/* compute sum's */
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	sDptr = sD->vec--;
	Nptr = N--;
	i = 1;
	while(i < D->n) {
		++sDptr;
		sDvec = sD->vec;
		++Nptr;
		Nvec = N;
		j = ++i;
		while(--j) {
			dist = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= dist) {
				*sDptr += dist;
				*++sDvec += dist;
				++*Nptr;
				++*++Nvec;
			} else {
				++sDvec;
				++Nvec;
			}
		}
	}
	++sD->vec;
	++N;
	
	return N;
}

long unsigned initQ(Matrix *D, Vector *sD, int *N) {
	
	int i, j, mi, mj, N_i, *Ni, *Nj;
	long unsigned pos;
	double min, q, sD_i, *sDi, *sDj, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* init */
	mi = 0;
	mj = 0;
	min = 1;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	sDi = sD->vec;
	Ni = N;
	i = 0;
	while(++i < D->n) {
		N_i = *++Ni;
		Nj = N - 1;
		sD_i = *++sDi;
		sDj = sD->vec - 1;
		j = -1;
		while(++j < i) {
			q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= q) {
				q = ((N_i + *++Nj - 4) >> 1) * q - sD_i - *++sDj;
				if(q <= min) {
					min = q;
					mi = i;
					mj = j;
				}
			} else {
				++Nj;
				++sDj;
			}
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

double initQchunk(Matrix *D, double *sD, int *N, double min, int *mi, int *mj, int i, int j) {
	
	const int chunk = 65536;
	int end, N_i;
	double sD_i, Q, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/* init */
	end = (j + chunk < i) ? (j + chunk) : i;
	N_i = N[i];
	N += --j;
	sD_i = sD[i];
	sD += j;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	
	/* get chunk */
	while(++j < end) {
		Q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
		if(0 <= Q) {
			Q = ((N_i + *++N - 4) >> 1) * Q - sD_i - *++sD;
			if(Q <= min) {
				min = Q;
				*mi = i;
				*mj = j;
			}
		} else {
			++N;
			++sD;
		}
	}
	
	return min;
}

long unsigned initQ_MN(Matrix *D, Vector *sD, int *N) {
	
	int i, j, mi, mj, N_i, *Ni, *Nj;
	long unsigned pos;
	double max, Q, sD_i, *sDi, *sDj, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* init */
	mi = 0;
	mj = 0;
	max = -DBL_MAX;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	sDi = sD->vec;
	Ni = N;
	i = 0;
	while(++i < D->n) {
		N_i = *++Ni;
		Nj = N - 1;
		sD_i = *++sDi;
		sDj = sD->vec - 1;
		j = -1;
		while(++j < i) {
			Q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= Q) {
				Q = ((N_i + *++Nj - 4) >> 1) * Q - sD_i - *++sDj;
				if(max <= Q) {
					max = Q;
					mi = i;
					mj = j;
				}
			} else {
				++Nj;
				++sDj;
			}
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

double initQ_MNchunk(Matrix *D, double *sD, int *N, double max, int *mi, int *mj, int i, int j) {
	
	const int chunk = 65536;
	int end, N_i;
	double sD_i, Q, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/* init */
	end = (j + chunk < i) ? (j + chunk) : i;
	N_i = N[i];
	N += --j;
	sD_i = sD[i];
	sD += j;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	
	/* get chunk */
	while(++j < end) {
		Q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
		if(0 <= Q) {
			Q = ((N_i + *++N - 4) >> 1) * Q - sD_i - *++sD;
			if(max <= Q) {
				max = Q;
				*mi = i;
				*mj = j;
			}
		} else {
			++N;
			++sD;
		}
	}
	
	return max;
}

void * initQ_thread(void *arg) {
	
	const int chunk = 65536;
	static volatile int thread_begin = 1, thread_wait, Lock = 0;
	static int thread_num = 1, next_i, next_j, Mi, Mj;
	static double Min;
	volatile int *lock = &Lock;
	NJthread *thread = arg;
	int i, j, mi, mj, *N;
	double min, *sDvec;
	Matrix *D;
	Vector *sD;
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* init */
	if(!thread) {
		wait_atomic((thread_begin != thread_num));
		thread_num = 0;
		thread_begin = 0;
		return NULL;
	}
	D = thread->D;
	sD = thread->sD;
	N = thread->N;
	if(thread->min) {
		lock(lock);
		wait_atomic((thread_begin != thread_num));
		/* init start */
		Mi = 0;
		Mj = 0;
		Min = (initQchunkPtr == &initQchunk) ? 1 : -DBL_MAX;
		next_i = 1;
		next_j = 0;
		thread_num = thread->num;
		thread_wait = thread_num;
		thread_begin = 0;
		unlock(lock);
	}
	do {
		wait_atomic(thread_begin);
		/* terminate */
		if(!thread_num) {
			return NULL;
		}

		/* init */
		i = 0;
		mi = 0;
		mj = 0;
		min = (initQchunkPtr == &initQchunk) ? 1 : -DBL_MAX;
		sDvec = sD->vec;

		/* fill in chunks */
		while(i < D->n) {
			/* get next chunk */
			lockTime(lock, 2);
			i = next_i;
			j = next_j;
			if(i <= (next_j += chunk)) {
				++next_i;
				next_j = 0;
			}
			unlock(lock);

			/* init chunk */
			if(i < D->n) {
				min = initQchunkPtr(D, sDvec, N, min, &mi, &mj, i, j);
			}
		}
		/* check min */
		lockTime(lock, 2);
		if(initQchunkPtr == &initQ_MNchunk) {
			if(Min < min) {
				Min = min;
				Mi = mi;
				Mj = mj;
			} else if(Min == min && (Mi < mi || (Mi == mi && Mj < mj))) {
				/* ensure reproducibility */
				Mi = mi;
				Mj = mj;
			}
		} else if(min < Min) {
			Min = min;
			Mi = mi;
			Mj = mj;
		} else if(min == Min && (Mi < mi || (Mi == mi && Mj < mj))) {
			/* ensure reproducibility */
			Mi = mi;
			Mj = mj;
		}
		--thread_wait;
		unlock(lock);
		wait_atomic(thread_wait);
		__sync_add_and_fetch(&thread_begin, 1);
	} while(!thread->min);
	thread->min = Min;
	thread->mi = Mi;
	thread->mj = Mj;

	return NULL;
}

long unsigned minD(Matrix *D, Vector *sD, int *N) {

	int i, j, mi, mj;
	long unsigned pos;
	double min, dist, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;

	mi = 0;
	mj = 0;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	min = (minDchunkPtr == &minDchunk) ? DBL_MAX : -DBL_MAX;
	for(i = 1; i < D->n; ++i) {
		for(j = 0; j < i; ++j) {
			dist = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= dist && dist <= min) {
				min = dist;
				mi = i;
				mj = j;
			}
		}
	}

	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;

	return pos;
}

double minDchunk(Matrix *D, int *N, double min, int *mi, int *mj, int i, int j) {

	const int chunk = 65536;
	int end;
	double dist, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;

	/* init */
	end = (j + chunk < i) ? (j + chunk) : i;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}

	/* get chunk */
	while(++j < end) {
		dist = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
		if(0 <= dist && dist <= min) {
			min = dist;
			*mi = i;
			*mj = j;
		}
	}

	return min;
}

long unsigned maxD(Matrix *D, Vector *sD, int *N) {

	int i, j, mi, mj;
	long unsigned pos;
	double max, dist, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;

	mi = 0;
	mj = 0;
	max = -DBL_MAX;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	for(i = 1; i < D->n; ++i) {
		for(j = 0; j < i; ++j) {
			dist = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= dist && max <= dist) {
				max = dist;
				mi = i;
				mj = j;
			}
		}
	}

	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;

	return pos;
}

double maxDchunk(Matrix *D, int *N, double max, int *mi, int *mj, int i, int j) {

	const int chunk = 65536;
	int end;
	double dist, *Dptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;

	/* init */
	end = (j + chunk < i) ? (j + chunk) : i;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}

	/* get chunk */
	while(++j < end) {
		dist = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
		if(0 <= dist && max <= dist) {
			max = dist;
			*mi = i;
			*mj = j;
		}
	}

	return max;
}

void * minD_thread(void *arg) {

	const int chunk = 65536;
	static volatile int thread_begin = 1, thread_wait, Lock = 0;
	static int thread_num = 1, next_i, next_j, Mi, Mj;
	static double Min;
	volatile int *lock = &Lock;
	NJthread *thread = arg;
	int i, j, mi, mj, *N;
	double min;
	Matrix *D;

	/* init */
	if(!thread) {
		wait_atomic((thread_begin != thread_num));
		thread_num = 0;
		thread_begin = 0;
		return NULL;
	}
	D = thread->D;
	N = thread->N;
	if(thread->min) {
		lock(lock);
		/* init start */
		Mi = 0;
		Mj = 0;
		Min = (minDchunkPtr == &minDchunk) ? DBL_MAX : -2;
		next_i = 1;
		next_j = 0;
		wait_atomic((thread_begin != thread_num));
		thread_num = thread->num;
		thread_wait = thread_num;
		thread_begin = 0;
		unlock(lock);
	}

	do {
		wait_atomic(thread_begin);
		/* terminate */
		if(!thread_num) {
			return NULL;
		}

		/* init */
		i = 0;
		mi = 0;
		mj = 0;
		min = (minDchunkPtr == &minDchunk) ? DBL_MAX : -2;

		/* fill in chunks */
		while(i < D->n) {
			/* get next chunk */
			lockTime(lock, 2);
			i = next_i;
			j = next_j;
			if(i <= (next_j += chunk)) {
				++next_i;
				next_j = 0;
			}
			unlock(lock);

			/* init chunk */
			if(i < D->n) {
				min = minDchunkPtr(D, N, min, &mi, &mj, i, j);
			}
		}

		/* check min */
		lockTime(lock, 2);
		if(minDchunkPtr == &maxDchunk) {
			if(Min < min) {
				Min = min;
				Mi = mi;
				Mj = mj;
			} else if(Min == min && (Mi < mi || (Mi == mi && Mj < mj))) {
				/* ensure reproducibility */
				Mi = mi;
				Mj = mj;
			}
		} else if(min < Min) {
			Min = min;
			Mi = mi;
			Mj = mj;
		} else if(min == Min && (Mi < mi || (Mi == mi && Mj < mj))) {
			/* ensure reproducibility */
			Mi = mi;
			Mj = mj;
		}
		--thread_wait;
		unlock(lock);
		wait_atomic(thread_wait);
		__sync_add_and_fetch(&thread_begin, 1);
	} while(!thread->min);

	thread->min = Min;
	thread->mi = Mi;
	thread->mj = Mj;

	return NULL;
}

long unsigned minPair(Matrix *Q) {

	int i, j, mi, mj;
	long unsigned pos;
	double min, q, *Qptr;
	float *Qfptr;
	short unsigned *Qsptr;
	unsigned char *Qbptr;

	mi = 0;
	mj = 0;
	min = 1;
	Qptr = 0;
	Qfptr = 0;
	Qsptr = 0;
	Qbptr = 0;
	if(Q->mat) {
		Qptr = *(Q->mat) - 1;
	} else if(Q->fmat) {
		Qfptr = *(Q->fmat) - 1;
	} else if(Q->smat) {
		Qsptr = *(Q->smat) - 1;
	} else {
		Qbptr = *(Q->bmat) - 1;
	}
	for(i = 1; i < Q->n; ++i) {
		for(j = 0; j < i; j++) {
			q = Qptr ? *++Qptr : Qfptr ? *++Qfptr : Qsptr ? uctod(*++Qsptr) : uctod(*++Qbptr);
			if(q <= min) {
				min = q;
				mi = i;
				mj = j;
			}
		}
	}

	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;

	return pos;
}

void updateD(Matrix *D, Vector *sD, int *N, int i, int j, double Li, double Lj) {

	int k, n, Dn, *Nptr;
	double dist, sd, D_ik, D_kj, D_ij, *D_i, *D_j, **Dmat, *sDvec;
	float *Df_i, *Df_j, **Dfmat;
	short unsigned *Ds_i, *Ds_j, **Dsmat;
	unsigned char *Db_i, *Db_j, **Dbmat;

	/*
	D(i,j) = D(j,i) <-> ltd -> j < i
	D(k,m) = (D(i,k) + D(k,j) - D(i,j))/2
	! D(i,k) -> D(k,m) = D(k,j) - L(j)
	! D(k,j) -> D(k,m) = D(i,k) - L(i)
	*/

	/* init */
	Dmat = 0;
	D_i = 0;
	D_j = 0;
	Dfmat = 0;
	Df_i = 0;
	Df_j = 0;
	Dsmat = 0;
	Ds_i = 0;
	Ds_j = 0;
	Dbmat = 0;
	Db_i = 0;
	Db_j = 0;
	if(D->mat) {
		Dmat = D->mat;
		D_i = Dmat[i] - 1;
		D_j = Dmat[j] - 1;
		D_ij = Dmat[i][j];
	} else if(D->fmat) {
		Dfmat = D->fmat;
		Df_i = Dfmat[i] - 1;
		Df_j = Dfmat[j] - 1;
		D_ij = Dfmat[i][j];
	} else if(D->smat) {
		Dsmat = D->smat;
		Ds_i = Dsmat[i] - 1;
		Ds_j = Dsmat[j] - 1;
		D_ij = uctod(Dsmat[i][j]);
	} else {
		Dbmat = D->bmat;
		Db_i = Dbmat[i] - 1;
		Db_j = Dbmat[j] - 1;
		D_ij = uctod(Dbmat[i][j]);
	}
	sDvec = sD->vec - 1;
	Nptr = N - 1;

	/* prepare upadate on N and sD */
	n = 1;
	sd = 0;

	/* update (first) row */
	k = j + 1;
	while(--k) {
		D_ik = D_i ? *++D_i : Df_i ? *++Df_i : Ds_i ? uctod(*++Ds_i) : uctod(*++Db_i);
		D_kj = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
		if(0 <= D_ik && 0 <= D_kj) {
			dist = (D_ik + D_kj - D_ij) / 2;
			dist = dist < 0 ? 0 : dist; /* hnj approx-error */
			if(D_j) {
				*D_j = dist;
			} else if(Df_j) {
				*Df_j = dist;
			} else if(Ds_j) {
				*Ds_j = dtouc(dist, 0.25);
			} else {
				*Db_j = dtouc(dist, 0.25);
			}
			/* update N and sD */
			*++sDvec -= (D_ik + D_kj - dist);
			sd += dist;
			--*++Nptr;
			++n;
		} else if(0 <= D_ik) {
			dist = D_ik - Li;
			if(D_j) {
				*D_j = dist;
			} else if(Df_j) {
				*Df_j = dist;
			} else if(Ds_j) {
				*Ds_j = dtouc(dist, 0);
			} else {
				*Db_j = dtouc(dist, 0);
			}
			/* update N and sD, sD(k) is new */
			*++sDvec -= Li;
			++Nptr;
			sd += dist;
			++n;
		} else if(0 <= D_kj) {
			if(D_j) {
				dist = (*D_j -= Lj);
			} else if(Df_j) {
				dist = (*Df_j -= Lj);
			} else if(Ds_j) {
				dist = (*Ds_j -= dtouc(Lj, 0));
				dist = uctod(dist);
			} else {
				dist = (*Db_j -= dtouc(Lj, 0));
				dist = uctod(dist);
			}
			/* update N and sD */
			*++sDvec += (dist - D_kj);
			--*++Nptr;
			sd += dist;
			++n;
		}
	}

	/* update jth column */
	if(D_j) {
		D_j = Dmat[j];
	} else if(Df_j) {
		Df_j = Dfmat[j];
	} else if(Ds_j) {
		Ds_j = Dsmat[j];
	} else {
		Db_j = Dbmat[j];
	}
	++sDvec;
	++Nptr;
	Dn = i;
	k = j;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++sDvec;
			++Nptr;
		}
		while(++k < Dn) {
			if(Dmat) {
				D_ik = k < i ? Dmat[i][k] : Dmat[k][i];
				D_kj = Dmat[k][j];
			} else if(Dfmat) {
				D_ik = k < i ? Dfmat[i][k] : Dfmat[k][i];
				D_kj = Dfmat[k][j];
			} else if(Dsmat) {
				D_ik = k < i ? Dsmat[i][k] : Dsmat[k][i];
				D_ik = uctod(D_ik);
				D_kj = uctod(Dsmat[k][j]);
			} else {
				D_ik = k < i ? Dbmat[i][k] : Dbmat[k][i];
				D_ik = uctod(D_ik);
				D_kj = uctod(Dbmat[k][j]);
			}
			if(0 <= D_ik && 0 <= D_kj) {
				dist = (D_kj + D_ik - D_ij) / 2;
				dist = dist < 0 ? 0 : dist; /* hnj approx-error */
				if(Dmat) {
					Dmat[k][j] = dist;
				} else if(Dfmat) {
					Dfmat[k][j] = dist;
				} else if(Dsmat) {
					Dsmat[k][j] = dtouc(dist, 0.25);
				} else {
					Dbmat[k][j] = dtouc(dist, 0.25);
				}
				/* update N and sD */
				*++sDvec -= (D_ik + D_kj - dist);
				--*++Nptr;
				sd += dist;
				++n;
			} else if(0 <= D_ik) {
				dist = D_ik - Li;
				if(Dmat) {
					Dmat[k][j] = dist;
				} else if(Dfmat) {
					Dfmat[k][j] = dist;
				} else if(Dsmat) {
					Dsmat[k][j] = dtouc(dist, 0);
				} else {
					Dbmat[k][j] = dtouc(dist, 0);
				}
				/* update N and sD, sD(k) is new */
				*++sDvec -= Li;
				++Nptr;
				sd += dist;
				++n;
			} else if(0 <= D_kj) {
				if(Dmat) {
					dist = (Dmat[k][j] -= Lj) - D_j[k];
				} else if(Dfmat) {
					dist = (Dfmat[k][j] -= Lj) - Df_j[k];
				} else if(Dsmat) {
					dist = (Dsmat[k][j] -= dtouc(Lj, 0)) - Ds_j[k];
					dist = uctod(dist);
				} else {
					dist = (Dbmat[k][j] -= dtouc(Lj, 0)) - Db_j[k];
					dist = uctod(dist);
				}
				/* update N and sD */
				*++sDvec += dist;
				--*++Nptr;
				sd += dist;
				++n;
			}
		}
	}

	/* update N and sD for row j*/
	N[j] = n;
	sD->vec[j] = sd;
}

void updateD_UPGMA(Matrix *D, Vector *sD, int *N, int i, int j, double Li, double Lj) {

	int k, n, Dn, *Nptr;
	double dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec;
	float *Df_i, *Df_j, **Dfmat;
	short unsigned *Ds_i, *Ds_j, **Dsmat;
	unsigned char *Db_i, *Db_j, **Dbmat;

	/*
	UPGMA
	*/

	/* init */
	Dmat = 0;
	D_i = 0;
	D_j = 0;
	Dfmat = 0;
	Df_i = 0;
	Df_j = 0;
	Dsmat = 0;
	Ds_i = 0;
	Ds_j = 0;
	Dbmat = 0;
	Db_i = 0;
	Db_j = 0;
	if(D->mat) {
		Dmat = D->mat;
		D_i = Dmat[i] - 1;
		D_j = Dmat[j] - 1;
	} else if(D->fmat) {
		Dfmat = D->fmat;
		Df_i = Dfmat[i] - 1;
		Df_j = Dfmat[j] - 1;
	} else if(D->smat) {
		Dsmat = D->smat;
		Ds_i = Dsmat[i] - 1;
		Ds_j = Dsmat[j] - 1;
	} else {
		Dbmat = D->bmat;
		Db_i = Dbmat[i] - 1;
		Db_j = Dbmat[j] - 1;
	}
	sDvec = sD->vec - 1;
	Nptr = N - 1;

	/* prepare upadate on N and sD */
	n = 1;
	sd = 0;

	/* update (first) row */
	k = j + 1;
	while(--k) {
		D_ik = D_i ? *++D_i : Df_i ? *++Df_i : Ds_i ? uctod(*++Ds_i) : uctod(*++Db_i);
		D_kj = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
		if(0 <= D_ik && 0 <= D_kj) {
			dist = (D_ik + D_kj) / 2;
			if(D_j) {
				*D_j = dist;
			} else if(Df_j) {
				*Df_j = dist;
			} else if(Ds_j) {
				*Ds_j = dtouc(dist, 0.25);
			} else {
				*Db_j = dtouc(dist, 0.25);
			}
			/* update N and sD */
			*++sDvec -= (D_ik + D_kj - dist);
			sd += dist;
			--*++Nptr;
			++n;
		} else if(0 <= D_ik) {
			if(D_j) {
				*D_j = D_ik;
			} else if(Df_j) {
				*Df_j = D_ik;
			} else if(Ds_j) {
				*Ds_j = dtouc(D_ik, 0);
			} else {
				*Db_j = dtouc(D_ik, 0);
			}
			/* update N and sD, sD(k) is new */
			++sDvec;
			++Nptr;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			/* update N and sD */
			++sDvec;
			--*++Nptr;
			sd += D_kj;
			++n;
		}
	}

	/* update jth column */
	if(D_j) {
		D_j = Dmat[j];
	} else if(Df_j) {
		Df_j = Dfmat[j];
	} else if(Ds_j) {
		Ds_j = Dsmat[j];
	} else {
		Db_j = Dbmat[j];
	}
	++sDvec;
	++Nptr;
	Dn = i;
	k = j;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++sDvec;
			++Nptr;
		}
		while(++k < Dn) {
			if(Dmat) {
				D_ik = k < i ? Dmat[i][k] : Dmat[k][i];
				D_kj = Dmat[k][j];
			} else if(Dfmat) {
				D_ik = k < i ? Dfmat[i][k] : Dfmat[k][i];
				D_kj = Dfmat[k][j];
			} else if(Dsmat) {
				D_ik = k < i ? Dsmat[i][k] : Dsmat[k][i];
				D_ik = uctod(D_ik);
				D_kj = uctod(Dsmat[k][j]);
			} else {
				D_ik = k < i ? Dbmat[i][k] : Dbmat[k][i];
				D_ik = uctod(D_ik);
				D_kj = uctod(Dbmat[k][j]);
			}
			if(0 <= D_ik && 0 <= D_kj) {
				dist = (D_kj + D_ik) / 2;
				if(Dmat) {
					Dmat[k][j] = dist;
				} else if(Dfmat) {
					Dfmat[k][j] = dist;
				} else if(Dsmat) {
					Dsmat[k][j] = dtouc(dist, 0.25);
				} else {
					Dbmat[k][j] = dtouc(dist, 0.25);
				}
				/* update N and sD */
				*++sDvec -= (D_ik + D_kj - dist);
				--*++Nptr;
				sd += dist;
				++n;
			} else if(0 <= D_ik) {
				if(Dmat) {
					Dmat[k][j] = D_ik;
				} else if(Dfmat) {
					Dfmat[k][j] = D_ik;
				} else if(Dsmat) {
					Dsmat[k][j] = dtouc(D_ik, 0);
				} else {
					Dbmat[k][j] = dtouc(D_ik, 0);
				}
				/* update N and sD, sD(k) is new */
				++sDvec;
				++Nptr;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				/* update N and sD */
				++sDvec;
				--*++Nptr;
				sd += D_kj;
				++n;
			}
		}
	}

	/* update N and sD for row j*/
	N[j] = n;
	sD->vec[j] = sd;
}

void updateD_FF(Matrix *D, Vector *sD, int *N, int i, int j, double Li, double Lj) {

	int k, n, Dn, *Nptr;
	double dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec;
	float *Df_i, *Df_j, **Dfmat;
	short unsigned *Ds_i, *Ds_j, **Dsmat;
	unsigned char *Db_i, *Db_j, **Dbmat;

	/*
	Furthest First
	*/

	/* init */
	Dmat = 0;
	D_i = 0;
	D_j = 0;
	Dfmat = 0;
	Df_i = 0;
	Df_j = 0;
	Dsmat = 0;
	Ds_i = 0;
	Ds_j = 0;
	Dbmat = 0;
	Db_i = 0;
	Db_j = 0;
	if(D->mat) {
		Dmat = D->mat;
		D_i = Dmat[i] - 1;
		D_j = Dmat[j] - 1;
	} else if(D->fmat) {
		Dfmat = D->fmat;
		Df_i = Dfmat[i] - 1;
		Df_j = Dfmat[j] - 1;
	} else if(D->smat) {
		Dsmat = D->smat;
		Ds_i = Dsmat[i] - 1;
		Ds_j = Dsmat[j] - 1;
	} else {
		Dbmat = D->bmat;
		Db_i = Dbmat[i] - 1;
		Db_j = Dbmat[j] - 1;
	}
	sDvec = sD->vec - 1;
	Nptr = N - 1;

	/* prepare upadate on N and sD */
	n = 1;
	sd = 0;

	/* update (first) row */
	k = j + 1;
	while(--k) {
		D_ik = D_i ? *++D_i : Df_i ? *++Df_i : Ds_i ? uctod(*++Ds_i) : uctod(*++Db_i);
		D_kj = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
		if(0 <= D_ik && 0 <= D_kj) {
			dist = D_ik < D_kj ? D_kj : D_ik;
			if(D_j) {
				*D_j = dist;
			} else if(Df_j) {
				*Df_j = dist;
			} else if(Ds_j) {
				*Ds_j = dtouc(dist, 0);
			} else {
				*Db_j = dtouc(dist, 0);
			}
			/* update N and sD */
			*++sDvec -= (D_ik + D_kj - dist);
			sd += dist;
			--*++Nptr;
			++n;
		} else if(0 <= D_ik) {
			if(D_j) {
				*D_j = D_ik;
			} else if(Df_j) {
				*Df_j = D_ik;
			} else if(Ds_j) {
				*Ds_j = dtouc(D_ik, 0);
			} else {
				*Db_j = dtouc(D_ik, 0);
			}
			/* update N and sD, sD(k) is new */
			++sDvec;
			++Nptr;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			/* update N and sD */
			++sDvec;
			--*++Nptr;
			sd += D_kj;
			++n;
		}
	}

	/* update jth column */
	++sDvec;
	++Nptr;
	Dn = i;
	k = j;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++sDvec;
			++Nptr;
		}
		while(++k < Dn) {
			if(Dmat) {
				D_ik = k < i ? Dmat[i][k] : Dmat[k][i];
				D_kj = Dmat[k][j];
			} else if(Dfmat) {
				D_ik = k < i ? Dfmat[i][k] : Dfmat[k][i];
				D_kj = Dfmat[k][j];
			} else if(Dsmat) {
				D_ik = k < i ? Dsmat[i][k] : Dsmat[k][i];
				D_ik = uctod(D_ik);
				D_kj = uctod(Dsmat[k][j]);
			} else {
				D_ik = k < i ? Dbmat[i][k] : Dbmat[k][i];
				D_ik = uctod(D_ik);
				D_kj = uctod(Dbmat[k][j]);
			}
			if(0 <= D_ik && 0 <= D_kj) {
				dist = D_ik < D_kj ? D_kj : D_ik;
				if(Dmat) {
					Dmat[k][j] = dist;
				} else if(Dfmat) {
					Dfmat[k][j] = dist;
				} else if(Dsmat) {
					Dsmat[k][j] = dtouc(dist, 0);
				} else {
					Dbmat[k][j] = dtouc(dist, 0);
				}
				/* update N and sD */
				*++sDvec -= (D_ik + D_kj - dist);
				sd += dist;
				--*++Nptr;
				++n;
			} else if(0 <= D_ik) {
				if(Dmat) {
					Dmat[k][j] = D_ik;
				} else if(Dfmat) {
					Dfmat[k][j] = D_ik;
				} else if(Dsmat) {
					Dsmat[k][j] = dtouc(D_ik, 0);
				} else {
					Dbmat[k][j] = dtouc(D_ik, 0);
				}
				/* update N and sD, sD(k) is new */
				++sDvec;
				++Nptr;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				/* update N and sD */
				++sDvec;
				--*++Nptr;
				sd += D_kj;
				++n;
			}
		}
	}

	/* update N and sD for row j*/
	N[j] = n;
	sD->vec[j] = sd;
}

void updateD_CF(Matrix *D, Vector *sD, int *N, int i, int j, double Li, double Lj) {

	int k, n, Dn, *Nptr;
	double dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec;
	float *Df_i, *Df_j, **Dfmat;
	short unsigned *Ds_i, *Ds_j, **Dsmat;
	unsigned char *Db_i, *Db_j, **Dbmat;

	/*
	Closest First
	*/

	/* init */
	Dmat = 0;
	D_i = 0;
	D_j = 0;
	Dfmat = 0;
	Df_i = 0;
	Df_j = 0;
	Dsmat = 0;
	Ds_i = 0;
	Ds_j = 0;
	Dbmat = 0;
	Db_i = 0;
	Db_j = 0;
	if(D->mat) {
		Dmat = D->mat;
		D_i = Dmat[i] - 1;
		D_j = Dmat[j] - 1;
	} else if(D->fmat) {
		Dfmat = D->fmat;
		Df_i = Dfmat[i] - 1;
		Df_j = Dfmat[j] - 1;
	} else if(D->smat) {
		Dsmat = D->smat;
		Ds_i = Dsmat[i] - 1;
		Ds_j = Dsmat[j] - 1;
	} else {
		Dbmat = D->bmat;
		Db_i = Dbmat[i] - 1;
		Db_j = Dbmat[j] - 1;
	}
	sDvec = sD->vec - 1;
	Nptr = N - 1;

	/* prepare upadate on N and sD */
	n = 1;
	sd = 0;

	/* update (first) row */
	k = j + 1;
	while(--k) {
		D_ik = D_i ? *++D_i : Df_i ? *++Df_i : Ds_i ? uctod(*++Ds_i) : uctod(*++Db_i);
		D_kj = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
		if(0 <= D_ik && 0 <= D_kj) {
			dist = D_ik < D_kj ? D_ik : D_kj;
			if(D_j) {
				*D_j = dist;
			} else if(Df_j) {
				*Df_j = dist;
			} else if(Ds_j) {
				*Ds_j = dtouc(dist, 0);
			} else {
				*Db_j = dtouc(dist, 0);
			}
			/* update N and sD */
			*++sDvec -= (D_ik + D_kj - dist);
			sd += dist;
			--*++Nptr;
			++n;
		} else if(0 <= D_ik) {
			if(D_j) {
				*D_j = D_ik;
			} else if(Df_j) {
				*Df_j = D_ik;
			} else if(Ds_j) {
				*Ds_j = dtouc(D_ik, 0);
			} else {
				*Db_j = dtouc(D_ik, 0);
			}
			/* update N and sD, sD(k) is new */
			++sDvec;
			++N;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			/* update N and sD */
			++sDvec;
			--*++Nptr;
			sd += D_kj;
			++n;
		}
	}

	/* update jth column */
	++sDvec;
	++Nptr;
	Dn = i;
	k = j;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++sDvec;
			++Nptr;
		}
		while(++k < Dn) {
			if(Dmat) {
				D_ik = k < i ? Dmat[i][k] : Dmat[k][i];
				D_kj = Dmat[k][j];
			} else if(Dfmat) {
				D_ik = k < i ? Dfmat[i][k] : Dfmat[k][i];
				D_kj = Dfmat[k][j];
			} else if(Dsmat) {
				D_ik = k < i ? Dsmat[i][k] : Dsmat[k][i];
				D_ik = uctod(D_ik);
				D_kj = uctod(Dsmat[k][j]);
			} else {
				D_ik = k < i ? Dbmat[i][k] : Dbmat[k][i];
				D_ik = uctod(D_ik);
				D_kj = uctod(Dbmat[k][j]);
			}
			if(0 <= D_ik && 0 <= D_kj) {
				dist = D_ik < D_kj ? D_ik : D_kj;
				dist = dist < 0 ? 0 : dist;
				if(Dmat) {
					Dmat[k][j] = dist;
				} else if(Dfmat) {
					Dfmat[k][j] = dist;
				} else if(Dsmat) {
					Dsmat[k][j] = dtouc(dist, 0);
				} else {
					Dbmat[k][j] = dtouc(dist, 0);
				}
				/* update N and sD */
				*++sDvec -= (D_ik + D_kj - dist);
				--*++Nptr;
				sd += dist;
				++n;
			} else if(0 <= D_ik) {
				if(Dmat) {
					Dmat[k][j] = D_ik;
				} else if(Dfmat) {
					Dfmat[k][j] = D_ik;
				} else if(Dsmat) {
					Dsmat[k][j] = dtouc(D_ik, 0);
				} else {
					Dbmat[k][j] = dtouc(D_ik, 0);
				}
				/* update N and sD, sD(k) is new */
				++sDvec;
				++N;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				/* update N and sD */
				++sDvec;
				--*++Nptr;
				sd += D_kj;
				++n;
			}
		}
	}

	/* update N and sD for row j*/
	N[j] = n;
	sD->vec[j] = sd;
}

int * nj(Matrix *D, Vector *sD, int *N, Qseqs **names) {

	int i, j, mask, shift;
	long unsigned pair;
	double Li, Lj;
	Qseqs *tmp;

	/* init */
	N = initSummaD(sD, D, N);
	mask = UINT_MAX;
	shift = 8 * sizeof(unsigned);

	/* get pairs */
	while(D->n != 2 && (pair = minDist(D, sD, N))) {
		j = pair & mask;
		i = pair >> shift;

		/* get limbs */
		limbLengthPtr(&Li, &Lj, i, j, sD, N, (D->mat ? D->mat[i][j] : D->fmat ? D->fmat[i][j] : D->smat ? uctod(D->smat[i][j]) : uctod(D->bmat[i][j])));

		/* form leaf */
		formNode(names[j], names[i], Lj, Li);

		/* update D and vectors */
		updateDptr(D, sD, N, i, j, Li, Lj);

		/* rearrange */
		ltdMatrix_popArrange(D, i);
		sD->vec[i] = sD->vec[--sD->n];
		N[i] = N[D->n];
		exchange(names[i], names[D->n], tmp);
		ltdMatrixShrink(D, D->size - 1);
	}

	if(D->n == 2) {
		formLastNodePtr(*names, names[1], (D->mat ? **(D->mat) : D->fmat ? **(D->fmat) : D->smat ? uctod(**(D->smat)) : uctod(**(D->bmat))));
	} else {
		/* form remaining nodes with undefined distance */
		while(D->n != 1) {
			/* form leaf */
			formLastNodePtr(*names, names[--D->n], -1.0);
		}
	}

	/* verify valid newick format */
	if(*((*names)->seq) != '(') {
		byteshift((*names)->seq, (*names)->len, '(');
	}

	return N;
}

int * nj_thread(Matrix *D, Vector *sD, int *N, Qseqs **names, int thread_num) {

	int i, j;
	double Li, Lj;
	NJthread *threads, *thread;
	Qseqs *tmp;

	if(thread_num == 1) {
		return nj(D, sD, N, names);
	}

	/* init */
	N = initSummaD(sD, D, N);
	j = 0;

	/* start threads */
	threads = 0;
	thread = threads;
	i = thread_num;
	while(i) {
		thread = smalloc(sizeof(NJthread));
		thread->D = D;
		thread->sD = sD;
		thread->N = N;
		thread->min = 0;
		thread->mi = 0;
		thread->mj = 0;
		thread->num = thread_num;
		thread->next = threads;
		threads = thread;
		if(--i && (errno = pthread_create(&thread->id, NULL, minDist_thread, thread))) {
			fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
			fprintf(stderr, "Will continue with %d threads.\n", thread_num - i);
			i = 0;
		}
	}
	thread->id = 0;
	thread->min = 1;

	/* get pairs */
	while(D->n != 2 && (minDist_thread(thread) || ((i = thread->mi) | (j = thread->mj)))) {
		/* get limbs */
		limbLengthPtr(&Li, &Lj, i, j, sD, N, (D->mat ? D->mat[i][j] : D->fmat ? D->fmat[i][j] : D->smat ? uctod(D->smat[i][j]) : uctod(D->bmat[i][j])));

		/* form leaf */
		formNode(names[j], names[i], Lj, Li);

		/* update D and vectors */
		updateDptr(D, sD, N, i, j, Li, Lj);

		/* rearrange */
		ltdMatrix_popArrange(D, i);
		sD->vec[i] = sD->vec[--sD->n];
		N[i] = N[D->n];
		exchange(names[i], names[D->n], tmp);
		ltdMatrixShrink(D, D->size - 1);
		threads->min = 1;
	}

	if(D->n == 2) {
		formLastNodePtr(*names, names[1], (D->mat ? **(D->mat) : D->fmat ? **(D->fmat) : D->smat ? uctod(**(D->smat)) : uctod(**(D->bmat))));
	} else {
		/* form remaining nodes with undefined distance */
		while(D->n != 1) {
			/* form leaf */
			formLastNodePtr(*names, names[--D->n], -1.0);
		}
	}

	/* verify valid newick format */
	if(*((*names)->seq) != '(') {
		byteshift((*names)->seq, (*names)->len, '(');
	}

	/* collect threads */
	minDist_thread(0);
	thread = threads->next;
	free(threads);
	while(thread) {
		threads = thread->next;
		if((errno = pthread_join(thread->id, NULL))) {
			ERROR();
		}
		free(thread);
		thread = threads;
	}
    if (minDist_thread != NULL) {
        minDist_thread = NULL;
    }
	return N;
}
