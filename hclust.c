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
#include "dnjheaders/hclust.h"
#include "dnjheaders/matrix.h"
#include "dnjheaders/nj.h"
#include "dnjheaders/nwck.h"
#include "dnjheaders/qseqs.h"
#include "dnjheaders/pherror.h"
#include "dnjheaders/str.h"
#include "dnjheaders/vector.h"

int * (*initDsDQN)(Matrix *, Vector *, Vector *, int *) = &initHNJ;
long unsigned (*pairQ)(Vector *, int *) = &minQ;
int (*updateDsDQNPtr)(Matrix *, Vector *, Vector *, int *, int*, int, int, double, double) = &updateHNJ;
int (*popArrangePtr)(Matrix *, Vector *, Vector *, int*, int*, int) = &HNJ_popArrange;

int * initAlloc(Matrix *D, Vector *sD, Vector *Q, int *N) {
	
	/* alloc */
	sD->n = D->n;
	Q->n = D->n;
	if(sD->size < sD->n) {
		vector_realloc(sD, sD->n);
		vector_realloc(Q, Q->n);
		if(!(N = realloc(N, 2 * D->n * sizeof(int)))) {
				ERROR();
		}
	}
	
	return N;
}

int * initHNJ(Matrix *D, Vector *sD, Vector *Q, int *N) {
	
	int i, j, pos, N_i, *Ni, *Nj, *P, *Ptr;
	double min, minD, d, q, sD_i, *sDi, *sDj, *Dptr, *Qptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/* alloc */
	N = initAlloc(D, sD, Q, N);
	P = N + D->n;
	
	/* init sD and N */
	N = initSummaD(sD, D, N);
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* init */
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
	sDi = --sD->vec;
	Qptr = Q->vec - 1;
	Ni = --N;
	Ptr = P - 1;
	i = -1;
	while(++i < D->n) {
		sD_i = *++sDi;
		sDj = sD->vec;
		N_i = *++Ni;
		Nj = N;
		min = DBL_MAX;
		minD = DBL_MAX;
		pos = 0;
		j = -1;
		while(++j < i) {
			d = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= d) {
				q = ((N_i + *++Nj - 4) >> 1) * d - sD_i - *++sDj;
				if(q <= min) {
					if(q < min || d <= minD) {
						min = q;
						minD = d;
						pos = j;
					}
				}
			} else {
				++Nj;
				++sDj;
			}
		}
		/* store best match */
		*++Qptr = min;
		*++Ptr = pos;
	}
	++sD->vec;
	++N;
	
	return N;
}

int * initHMN(Matrix *D, Vector *sD, Vector *Q, int *N) {
	
	int i, j, pos, N_i, *Ni, *Nj, *P, *Ptr;
	double max, q, sD_i, *sDi, *sDj, *Dptr, *Qptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/* alloc */
	N = initAlloc(D, sD, Q, N);
	P = N + D->n;
	
	/* init sD and N */
	N = initSummaD(sD, D, N);
	
	/*
	Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	sD(i) = summa(d(i,k));
	org:
	Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
	*/
	
	/* init */
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
	sDi = --sD->vec;
	Qptr = Q->vec - 1;
	Ni = --N;
	Ptr = P - 1;
	i = -1;
	while(++i < D->n) {
		N_i = *++Ni;
		Nj = N;
		sD_i = *++sDi;
		sDj = sD->vec;
		max = -DBL_MAX;
		pos = 0;
		j = -1;
		while(++j < i) {
			q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= q) {
				q = q * ((N_i + *++Nj - 4) >> 1) - sD_i - *++sDj;
				if(max <= q) {
					max = q;
					pos = j;
					
				}
			} else {
				++Nj;
				++sDj;
			}
		}
		/* store best match */
		*++Qptr = max;
		*++Ptr = pos;
	}
	++sD->vec;
	++N;
	
	return N;
}

int * initDmin(Matrix *D, Vector *sD, Vector *Q, int *N) {
	
	int i, j, pos, *Nptr, *Nvec, *P, *Ptr;
	double min, dist, *sDptr, *sDvec, *Dptr, *Qptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/* alloc */
	N = initAlloc(D, sD, Q, N);
	P = N + D->n;
	
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
	sDptr = --sD->vec;
	Qptr = Q->vec - 1;
	Nptr = --N;
	Ptr = P - 1;
	i = -1;
	while(++i < D->n) {
		++sDptr;
		sDvec = sD->vec;
		++Nptr;
		Nvec = N;
		min = DBL_MAX;
		pos = 0;
		j = -1;
		while(++j < i) {
			dist = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= dist) {
				if(dist <= min) {
					min = dist;
					pos = j;
				}
				*sDptr += dist;
				*++sDvec += dist;
				++*Nptr;
				++*++Nvec;
			} else {
				++sDvec;
				++Nvec;
			}
		}
		/* store best match */
		*++Qptr = min;
		*++Ptr = pos;
	}
	++sD->vec;
	++N;
	
	return N;
}

int * initDmax(Matrix *D, Vector *sD, Vector *Q, int *N) {
	
	int i, j, pos, *Nptr, *Nvec, *P, *Ptr;
	double max, dist, *sDptr, *sDvec, *Dptr, *Qptr;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/* alloc */
	N = initAlloc(D, sD, Q, N);
	P = N + D->n;
	
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
	sDptr = --sD->vec;
	Qptr = Q->vec - 1;
	Nptr = --N;
	Ptr = P - 1;
	i = -1;
	while(++i < D->n) {
		++sDptr;
		sDvec = sD->vec;
		++Nptr;
		Nvec = N;
		pos = 0;
		max = -DBL_MAX;
		j = -1;
		while(++j < i) {
			dist = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
			if(0 <= dist) {
				if(max <= dist) {
					max = dist;
					pos = j;
				}
				*sDptr += dist;
				*++sDvec += dist;
				++*Nptr;
				++*++Nvec;
			} else {
				++sDvec;
				++Nvec;
			}
		}
		/* store best match */
		*++Qptr = max;
		*++Ptr = pos;
	}
	++sD->vec;
	++N;
	
	return N;
}

long unsigned minQ(Vector *Q, int *P) {
	
	int i, mi, mj;
	long unsigned pos;
	double min, *Qptr;
	
	mi = 0;
	mj = 0;
	min = DBL_MAX;
	Qptr = Q->vec;
	i = 0;
	while(++i < Q->n) {
		if(*++Qptr <= min) {
			min = *Qptr;
			mi = i;
			mj = *++P;
		} else {
			++P;
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

long unsigned maxQ(Vector *Q, int *P) {
	
	int i, mi, mj;
	long unsigned pos;
	double max, *Qptr;
	
	mi = 0;
	mj = 0;
	max = -DBL_MAX;
	Qptr = Q->vec;
	i = 0;
	while(++i < Q->n) {
		if(max <= *++Qptr) {
			max = *Qptr;
			mi = i;
			mj = *++P;
		} else {
			++P;
		}
	}
	
	/* save pos */
	pos = 0;
	pos |= mi;
	pos <<= sizeof(unsigned) * 8;
	pos |= mj;
	
	return pos;
}

void updatePrevQ(Matrix *D, Vector *sD, Vector *Q, int *N, int *P) {
	
	int i, *Nptr;
	double **Dmat, *sDvec, *Qvec, dist;
	float **Dfmat;
	short unsigned **Dsmat;
	unsigned char **Dbmat;
	
	/* update previous Q */
	Dmat = 0;
	Dfmat = 0;
	Dsmat = 0;
	Dbmat = 0;
	if(D->mat) {
		Dmat = D->mat - 1;
	} else if(D->fmat) {
		Dfmat = D->fmat - 1;
	} else if(D->smat) {
		Dsmat = D->smat - 1;
	} else {
		Dbmat = D->bmat - 1;
	}
	sDvec = sD->vec - 1;
	Qvec = Q->vec - 1;
	Nptr = N - 1;
	--P;
	i = D->n;
	while(--i) {
		dist = Dmat ? (*++Dmat)[*++P] : Dfmat ? (*++Dfmat)[*++P] : Dsmat ? uctod((*++Dsmat)[*++P]) : uctod((*++Dbmat)[*++P]);
		if(0 <= dist) {
			*++Qvec = ((*++Nptr + N[*P] - 4) >> 1) * dist - *++sDvec - sD->vec[*P];
		} else {
			++Qvec;
			++sDvec;
			++Nptr;
		}
	}
}

int updateHNJ(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj) {
	
	int k, p, Dn, N_j;
	double q, min, sD_j, *D_j, **Dmat, *sDvec, *Qvec;
	float *Df_j, **Dfmat;
	short unsigned *Ds_j, **Dsmat;
	unsigned char *Db_j, **Dbmat;
	
	/* update D, sD and N */
	updateD(D, sD, N, i, j, Li, Lj);
	
	/* update previous Q */
	updatePrevQ(D, sD, Q, N, P);
	
	/* update on (first) row */
	Qvec = Q->vec + j;
	*Qvec = DBL_MAX;
	P += j;
	*P = 0;
	Dmat = 0;
	D_j = 0;
	Dfmat = 0;
	Df_j = 0;
	Dsmat = 0;
	Ds_j = 0;
	Dbmat = 0;
	Db_j = 0;
	if(D->mat) {
		Dmat = D->mat;
		D_j = Dmat[j] - 1;
	} else if(D->fmat) {
		Dfmat = D->fmat;
		Df_j = Dfmat[j] - 1;
	} else if(D->smat) {
		Dsmat = D->smat;
		Ds_j = Dsmat[j] - 1;
	} else {
		Dbmat = D->bmat;
		Db_j = Dbmat[j] - 1;
	}
	sDvec = sD->vec - 1;
	sD_j = sD->vec[j];
	N_j = N[j];
	--N;
	for(k = 0; k < j; ++k) {
		q = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
		if(0 <= q) {
			q = ((N_j + *++N - 4) >> 1) * q - sD_j - *++sDvec;
			if(q <= *Qvec) {
				*Qvec = q;
				*P = k;
			}
		} else {
			++sDvec;
			++N;
		}
	}
	
	/* save min */
	min = *Qvec;
	p = j;
	
	/* update on jth column */
	++sDvec;
	++N;
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++Qvec;
			++P;
			++sDvec;
			++N;
		}
		while(++k < Dn) {
			q = Dmat ? Dmat[k][j] : Dfmat ? Dfmat[k][j] : Dsmat ? uctod(Dsmat[k][j]) : uctod(Dbmat[k][j]);
			if(0 <= q) {
				q = ((N_j + *++N - 4) >> 1) * q - sD_j - *++sDvec;
				/* check if best join was i or j, or new node is better */
				++P;
				++Qvec;
				if(*P == i || *P == j) {
					*Qvec = q;
					*P = j;
					if(q <= min) {
						q = min;
						p = k;
					}
				} else if(q <= *Qvec) {
					*Qvec = q;
					if(*P < j) {
						*P = j;
					}
					if(q <= min) {
						q = min;
						p = k;
					}
				}
			} else {
				++sDvec;
				++Qvec;
				++N;
				++P;
			}
		}
	}
	
	return p;
}

int updateHMN(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj) {
	
	int k, p, Dn, N_j;
	double q, max, sD_j, *D_j, **Dmat, *sDvec, *Qvec;
	float *Df_j, **Dfmat;
	short unsigned *Ds_j, **Dsmat;
	unsigned char *Db_j, **Dbmat;
	
	/* update D, sD and N */
	updateD(D, sD, N, i, j, Li, Lj);
	
	/* update previous Q */
	updatePrevQ(D, sD, Q, N, P);
	
	/* update on (first) row */
	Qvec = Q->vec + j;
	P += j;
	*Qvec = -DBL_MAX;
	Dmat = 0;
	D_j = 0;
	Dfmat = 0;
	Df_j = 0;
	Dsmat = 0;
	Ds_j = 0;
	Dbmat = 0;
	Db_j = 0;
	if(D->mat) {
		Dmat = D->mat;
		D_j = Dmat[j] - 1;
	} else if(D->fmat) {
		Dfmat = D->fmat;
		Df_j = Dfmat[j] - 1;
	} else if(D->smat) {
		Dsmat = D->smat;
		Ds_j = Dsmat[j] - 1;
	} else {
		Dbmat = D->bmat;
		Db_j = Dbmat[j] - 1;
	}
	sDvec = sD->vec - 1;
	sD_j = sD->vec[j];
	N_j = N[j];
	--N;
	for(k = 0; k < j; ++k) {
		q = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
		if(0 <= q) {
			q = ((N_j + *++N - 4) >> 1) * q - sD_j - *++sDvec;
			if(*Qvec < q) {
				*Qvec = q;
				*P = k;
			}
		} else {
			++sDvec;
			++N;
		}
	}
	
	/* save max */
	max = *Qvec;
	p = j;
	
	/* update on jth column */
	++sDvec;
	++N;
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++Qvec;
			++P;
			++sDvec;
			++N;
		}
		while(++k < Dn) {
			q = Dmat ? Dmat[k][j] : Dfmat ? Dfmat[k][j] : Dsmat ? uctod(Dsmat[k][j]) : uctod(Dbmat[k][j]);
			if(0 <= q) {
				q = ((N_j + *++N - 4) >> 1) * q - sD_j - *++sDvec;
				/* check if best join was i or j, or new node is better */
				++P;
				if(*++Qvec <= q || *P == i || *P == j) {
					if(*Qvec < q || *P == i || *P < j) {
						*P = j;
					}
					*Qvec = q;
					if(max <= q) {
						q = max;
						p = k;
					}
				}
			} else {
				++sDvec;
				++Qvec;
				++N;
				++P;
			}
		}
	}
	
	return p;
}

int updateUPGMA(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj) {
	
	int k, n, p, Dn, *Nptr;
	double min, dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec, *Qvec;
	float *Df_i, *Df_j, **Dfmat;
	short unsigned *Ds_i, *Ds_j, **Dsmat;
	unsigned char *Db_i, *Db_j, **Dbmat;
	
	/*
	UPGMA
	*/
	
	/* init */
	Qvec = Q->vec + j;
	*Qvec = DBL_MAX;
	P += j;
	*P = 0;
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
	for(k = 0; k < j; ++k) {
		D_ik = D_i ? *++D_i : Df_i ? *++Df_i : Ds_i ? uctod(*++Ds_i) : uctod(*++Db_i);
		D_kj = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
		if(0 <= D_ik && 0 <= D_kj) {
			dist = (D_ik + D_kj) / 2;
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
			dist = D_ik;
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
			++sDvec;
			++Nptr;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			dist = D_kj;
			/* update N and sD */
			++sDvec;
			--*++Nptr;
			sd += D_kj;
			++n;
		} else {
			dist = -1;
		}
		
		/* update Q and P */
		if(0 <= dist && dist <= *Qvec) {
			*Qvec = dist;
			*P = k;
		}
	}
	
	/* save min */
	min = *Qvec;
	p = j;
	
	/* update jth column */
	++sDvec;
	++Nptr;
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++sDvec;
			++Qvec;
			++Nptr;
			++P;
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
				dist = D_ik;
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
				++sDvec;
				++Nptr;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				dist = D_kj;
				/* update N and sD */
				++sDvec;
				--*++Nptr;
				sd += D_kj;
				++n;
			} else {
				dist = -1;
			}
			
			/* update Q and P */
			++Qvec;
			++P;
			if(0 <= dist) {
				if(dist < *Qvec) {
					/* dist < *Qvec -> k < i */
					*Qvec = dist;
					*P = j;
					if(min <= dist) {
						min = dist;
						p = k;
					}
				} else if(*P == i || *P == j) {
					/* old dist will be <= new dist */
					/* new dist will be >= old if *P == i v j */
					if(dist == *Qvec) {
						*P = j; /* d(i,k) == d(j,k) */
						if(min <= dist) {
							min = dist;
							p = k;
						}
					} else {
						*P = -1; /* mark as bounded row */
					}
				}
			}
		}
	}
	
	/* update N and sD for row j*/
	N[j] = n;
	sD->vec[j] = sd;
	
	return p;
}

int updateFF(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj) {
	
	int k, n, p, Dn, *Nptr;
	double min, dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec, *Qvec;
	float *Df_i, *Df_j, **Dfmat;
	short unsigned *Ds_i, *Ds_j, **Dsmat;
	unsigned char *Db_i, *Db_j, **Dbmat;
	
	/*
	Furthest First
	*/
	
	/* init */
	Qvec = Q->vec + j;
	P += j;
	*Qvec = DBL_MAX;
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
	for(k = 0; k < j; ++k) {
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
			dist = D_ik;
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
			++sDvec;
			++Nptr;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			dist = D_kj;
			/* update N and sD */
			++sDvec;
			--*++Nptr;
			sd += D_kj;
			++n;
		} else {
			dist = -1;
		}
		
		/* update Q and P */
		if(dist < *Qvec) {
			*Qvec = dist;
			*P = k;
		}
	}
	
	/* save min */
	min = *Qvec;
	p = j;
	
	/* update jth column */
	++sDvec;
	++Nptr;
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++Qvec;
			++P;
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
				dist = D_ik;
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
				++sDvec;
				++Nptr;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				dist = D_kj;
				/* update N and sD */
				++sDvec;
				--*++Nptr;
				sd += D_kj;
				++n;
			} else {
				dist = -1;
			}
			
			/* update Q and P */
			++Qvec;
			++P;
			if(0 <= dist) {
				if(dist < *Qvec) {
					/* dist < *Qvec -> k < i */
					*Qvec = dist;
					*P = j;
					if(min <= dist) {
						min = dist;
						p = k;
					}
				} else if(*P == i || *P == j) {
					/* old dist will be <= new dist */
					/* new dist will be >= old if *P == i v j */
					if(dist == *Qvec) {
						*P = j; /* d(i,k) == d(j,k) */
						if(min <= dist) {
							min = dist;
							p = k;
						}
					} else {
						*P = -1; /* mark as bounded row */
					}
				}
			}
		}
	}
	
	/* update N and sD for row j*/
	N[j] = n;
	sD->vec[j] = sd;
	
	return p;
}

int updateCF(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj) {
	
	int k, p, n, Dn, *Nptr;
	double min, dist, sd, D_ik, D_kj, *D_i, *D_j, **Dmat, *sDvec, *Qvec;
	float *Df_i, *Df_j, **Dfmat;
	short unsigned *Ds_i, *Ds_j, **Dsmat;
	unsigned char *Db_i, *Db_j, **Dbmat;
	
	/*
	Closest First
	*/
	
	/* init */
	Qvec = Q->vec + j;
	*Qvec = DBL_MAX;
	P += j;
	*P = 0;
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
	for(k = 0; k < j; ++k) {
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
			dist = D_ik;
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
			++sDvec;
			++N;
			sd += D_ik;
			++n;
		} else if(0 <= D_kj) {
			dist = D_kj;
			/* update N and sD */
			++sDvec;
			--*++Nptr;
			sd += D_kj;
			++n;
		} else {
			dist = -1;
		}
		
		/* update Q and P */
		if(0 <= dist && dist <= *Qvec) {
			*Qvec = dist;
			*P = k;
		}
	}
	
	/* save min */
	min = *Qvec;
	p = j;
	
	/* update jth column */
	++sDvec;
	++Nptr;
	Dn = i;
	while(Dn != D->n) {
		if(k == Dn) {
			/* skip ith row */
			Dn = D->n;
			++sDvec;
			++Qvec;
			++Nptr;
			++P;
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
				dist = D_ik;
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
				++sDvec;
				++N;
				sd += D_ik;
				++n;
			} else if(0 <= D_kj) {
				dist = D_kj;
				/* update N and sD */
				++sDvec;
				--*++Nptr;
				sd += D_kj;
				++n;
			} else {
				dist = -1;
			}
			
			/* update Q and P */
			++Qvec;
			++P;
			if(0 <= dist && dist <= *Qvec) {
				if(dist < *Qvec || *P == i || *P == k || *P < j) {
					*Qvec = dist;
					*P = j;
					if(min <= dist) {
						min = dist;
						p = k;
					}
				}
			}
		}
	}
	
	/* update N and sD for row j*/
	N[j] = n;
	sD->vec[j] = sd;
	
	return p;
}

int HNJ_popArrange(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int pos) {
	
	int i, n, p, N_i, *Nptr;
	double q, min, sD_i, *dest, *src, **mat, *sDvec, *Qptr;
	float *fdest, *fsrc, **fmat;
	short unsigned *sdest, *ssrc, **smat;
	unsigned char *bdest, *bsrc, **bmat;
	
	n = --D->n;
	--sD->n;
	--Q->n;
	if(pos != n) {
		sD->vec[pos] = sD->vec[n];
		N[pos] = N[n];
		Qptr = Q->vec + pos;
		q = DBL_MAX;
		*Qptr = q;
		P += pos;
		*P = 0;
		
		mat = 0;
		fmat = 0;
		smat = 0;
		bmat = 0;
		dest = 0;
		fdest = 0;
		sdest = 0;
		bdest = 0;
		src = 0;
		fsrc = 0;
		ssrc = 0;
		bsrc = 0;
		if(D->mat) {
			mat = D->mat;
			/* row to be emptied */
			dest = mat[pos] - 1;
			/* row to be moved up */
			src = mat[n] - 1;
		} else if(D->fmat) {
			fmat = D->fmat;
			/* row to be emptied */
			fdest = fmat[pos] - 1;
			/* row to be moved up */
			fsrc = fmat[n] - 1;
		} else if(D->smat) {
			smat = D->smat;
			/* row to be emptied */
			sdest = smat[pos] - 1;
			/* row to be moved up */
			ssrc = smat[n] - 1;
		} else {
			bmat = D->bmat;
			/* row to be emptied */
			bdest = bmat[pos] - 1;
			/* row to be moved up */
			bsrc = bmat[n] - 1;
		}
		
		/* copy last row into "first" row */
		sDvec = sD->vec - 1;
		sD_i = sD->vec[pos];
		Nptr = N - 1;
		N_i = N[pos];
		i = -1;
		while(++i < pos) {
			q = dest ? (*++dest = *++src) : fdest ? (*++fdest = *++fsrc) : sdest ? uctod((*++sdest = *++ssrc)) : uctod((*++bdest = *++bsrc));
			if(0 <= q) {
				q = q * ((N_i + *++Nptr - 4) >> 1) - sD_i - *++sDvec;
				if(q <= *Qptr) {
					*Qptr = q;
					*P = i;
				}
			} else {
				++sDvec;
				++Nptr;
			}
		}
		
		/* save min */
		min = *Qptr;
		p = pos;
		
		/* skip pos, as that is removied */
		if(src) {
			++src;
		} else if(fsrc) {
			++fsrc;
		} else if(ssrc) {
			++ssrc;
		} else {
			++bsrc;
		}
		++sDvec;
		++Nptr;
		
		/* tilt remaining part of last row into column "pos" */
		while(++i < n) {
			q = mat ? (mat[i][pos] = *++src) : fmat ? (fmat[i][pos] = *++fsrc) : smat ? uctod((smat[i][pos] = *++ssrc)) : uctod((bmat[i][pos] = *++bsrc));
			if(0 <= q) {
				q = q * ((N_i + *++Nptr - 4) >> 1) - sD_i - *++sDvec;
				if(q <= *++Qptr) {
					if(*++P < pos || q < *Qptr) {
						*Qptr = q;
						*P = pos;
						if(q <= min) {
							min = q;
							p = i;
						}
					}
				} else {
					++P;
				}
			} else {
				++sDvec;
				++Qptr;
				++Nptr;
				++P;
			}
		}
	} else {
		p = 0;
	}
	
	return p;
}

int HMN_popArrange(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int pos) {
	
	int i, n, p, N_i, *Nptr;
	double q, max, sD_i, *dest, *src, **mat, *sDvec, *Qptr;
	float *fdest, *fsrc, **fmat;
	short unsigned *sdest, *ssrc, **smat;
	unsigned char *bdest, *bsrc, **bmat;
	
	n = --D->n;
	--sD->n;
	--Q->n;
	if(pos != n) {
		sD->vec[pos] = sD->vec[n];
		N[pos] = N[n];
		Qptr = Q->vec + pos;
		*Qptr = -DBL_MAX;
		P += pos;
		*P = 0;
		
		mat = 0;
		fmat = 0;
		smat = 0;
		bmat = 0;
		dest = 0;
		fdest = 0;
		sdest = 0;
		bdest = 0;
		src = 0;
		fsrc = 0;
		ssrc = 0;
		bsrc = 0;
		if(D->mat) {
			mat = D->mat;
			/* row to be emptied */
			dest = mat[pos] - 1;
			/* row to be moved up */
			src = mat[n] - 1;
		} else if(D->fmat) {
			fmat = D->fmat;
			/* row to be emptied */
			fdest = fmat[pos] - 1;
			/* row to be moved up */
			fsrc = fmat[n] - 1;
		} else if(D->smat) {
			smat = D->smat;
			/* row to be emptied */
			sdest = smat[pos] - 1;
			/* row to be moved up */
			ssrc = smat[n] - 1;
		} else {
			bmat = D->bmat;
			/* row to be emptied */
			bdest = bmat[pos] - 1;
			/* row to be moved up */
			bsrc = bmat[n] - 1;
		}
		
		/* copy last row into "first" row */
		sDvec = sD->vec - 1;
		sD_i = sD->vec[pos];
		Nptr = N - 1;
		N_i = N[pos];
		i = -1;
		while(++i < pos) {
			q = dest ? (*++dest = *++src) : fdest ? (*++fdest = *++fsrc) : sdest ? uctod((*++sdest = *++ssrc)) : uctod((*++bdest = *++bsrc));
			if(0 <= q) {
				q = q * ((N_i + *++Nptr - 4) >> 1) - sD_i - *++sDvec;
				if(*Qptr <= q) {
					*Qptr = q;
					*P = i;
				}
			} else {
				++sDvec;
				++Nptr;
			}
		}
		
		/* save max */
		max = *Qptr;
		p = pos;
		
		/* skip pos, as that is removied */
		if(src) {
			++src;
		} else if(fsrc) {
			++fsrc;
		} else if(ssrc) {
			++ssrc;
		} else {
			++bsrc;
		}
		++sDvec;
		++Nptr;
		
		/* tilt remaining part of last row into column "pos" */
		while(++i < n) {
			q = mat ? (mat[i][pos] = *++src) : fmat ? (fmat[i][pos] = *++fsrc) : smat ? uctod((smat[i][pos] = *++ssrc)) : uctod((bmat[i][pos] = *++bsrc));
			if(0 <= q) {
				q = q * ((N_i + *++Nptr - 4) >> 1) - sD_i - *++sDvec;
				if(*++Qptr <= q) {
					if(*++P < pos || *Qptr < q) {
						*Qptr = q;
						*P = pos;
						if(max <= q) {
							q = max;
							p = i;
						}
					}
				} else {
					++P;
				}
			} else {
				++sDvec;
				++Qptr;
				++Nptr;
				++P;
			}
		}
	} else {
		p = 0;
	}
	
	return p;
}

int UPGMA_popArrange(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int pos) {
	
	/* same as for CF */
	int i, n, p;
	double q, min, *dest, *src, **mat, *Qptr;
	float *fdest, *fsrc, **fmat;
	short unsigned *sdest, *ssrc, **smat;
	unsigned char *bdest, *bsrc, **bmat;
	
	n = --D->n;
	--sD->n;
	--Q->n;
	if(pos != n) {
		sD->vec[pos] = sD->vec[n];
		N[pos] = N[n];
		Qptr = Q->vec + pos;
		q = DBL_MAX;
		*Qptr = q;
		P += pos;
		*P = 0;
		
		mat = 0;
		fmat = 0;
		smat = 0;
		bmat = 0;
		dest = 0;
		fdest = 0;
		sdest = 0;
		bdest = 0;
		src = 0;
		fsrc = 0;
		ssrc = 0;
		bsrc = 0;
		if(D->mat) {
			mat = D->mat;
			/* row to be emptied */
			dest = mat[pos] - 1;
			/* row to be moved up */
			src = mat[n] - 1;
		} else if(D->fmat) {
			fmat = D->fmat;
			/* row to be emptied */
			fdest = fmat[pos] - 1;
			/* row to be moved up */
			fsrc = fmat[n] - 1;
		} else if(D->smat) {
			smat = D->smat;
			/* row to be emptied */
			sdest = smat[pos] - 1;
			/* row to be moved up */
			ssrc = smat[n] - 1;
		} else {
			bmat = D->bmat;
			/* row to be emptied */
			bdest = bmat[pos] - 1;
			/* row to be moved up */
			bsrc = bmat[n] - 1;
		}
		
		/* copy last row into "first" row */
		i = -1;
		while(++i < pos) {
			q = dest ? (*++dest = *++src) : fdest ? (*++fdest = *++fsrc) : sdest ? uctod((*++sdest = *++ssrc)) : uctod((*++bdest = *++bsrc));
			if(0 <= q && q <= *Qptr) {
				*Qptr = q;
				*P = i;
			}
		}
		
		/* save min */
		min = *Qptr;
		p = pos;
		
		/* skip pos, as that is removied */
		if(src) {
			++src;
		} else if(fsrc) {
			++fsrc;
		} else if(ssrc) {
			++ssrc;
		} else {
			++bsrc;
		}
		
		/* tilt remaining part of last row into column "pos" */
		while(++i < n) {
			q = mat ? (mat[i][pos] = *++src) : fmat ? (fmat[i][pos] = *++fsrc) : smat ? uctod((smat[i][pos] = *++ssrc)) : uctod((bmat[i][pos] = *++bsrc));
			if(0 <= q) {
				if(q <= *++Qptr) {
					if(*++P < pos || q < *Qptr) {
						*Qptr = q;
						*P = pos;
						if(q <= min) {
							min = q;
							p = i;
						}
					}
				} else {
					++P;
				}
			} else {
				++Qptr;
				++P;
			}
		}
	} else {
		p = 0;		
	}
	
	return p;
}

int * hclust(Matrix *D, Vector *sD, Vector *Q, int *N, Qseqs **names) {
	
	int i, j, mask, shift, *P;
	long unsigned pair;
	double Li, Lj;
	Qseqs *tmp;
	
	/* init */
	N = initDsDQN(D, sD, Q, N);
	P = N + D->n;
	mask = UINT_MAX;
	shift = 8 * sizeof(unsigned);
	
	/* get pairs */
	while(D->n != 2 && (pair = pairQ(Q, P))) {
		j = pair & mask;
		i = pair >> shift;
		
		/* get limbs */
		limbLengthPtr(&Li, &Lj, i, j, sD, N, (D->mat ? D->mat[i][j] : D->fmat ? D->fmat[i][j] : D->smat ? uctod(D->smat[i][j]) : uctod(D->bmat[i][j])));
		
		/* form leaf */
		formNode(names[j], names[i], Lj, Li);
		
		/* update D and vectors */
		updateDsDQNPtr(D, sD, Q, N, P, i, j, Li, Lj);
		
		/* rearrange */
		popArrangePtr(D, sD, Q, N, P, i);
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
