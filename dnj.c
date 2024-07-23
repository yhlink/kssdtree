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

#define _XOPEN_SOURCE 600

#include <float.h>
#include <limits.h>
#include <pthread.h>
#include <stdlib.h>
#include "dnjheaders/bytescale.h"
#include "dnjheaders/dnj.h"
#include "dnjheaders/hclust.h"
#include "dnjheaders/matrix.h"
#include "dnjheaders/nj.h"
#include "dnjheaders/nwck.h"
#include "dnjheaders/qseqs.h"
#include "dnjheaders/pherror.h"
#include "dnjheaders/str.h"
#include "dnjheaders/threader.h"
#include "dnjheaders/vector.h"

long unsigned (*Qpair)(Matrix *, Vector *, Vector *, int *, int *, int) = &minQpair;

long unsigned (*Qrow)(Matrix *, Vector *, Vector *, int *, int *, double *, int) = minQrow;

int (*Qbool)(double, double, long unsigned, long unsigned) = &minQbool;

int (*nextQrow)(int, double, double *) = &nextQminRow;

int (*qPos)(double *, int, int) = &minPos;

long unsigned minQpair(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i) {

    int j, mj, N_i, *Nj;
    long unsigned pos;
    double min, q, updateQ, sD_i, *sDj, *Dptr, *Qptr;
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
    pos = 0;
    min = DBL_MAX;
    if (i && min != Q->vec[i]) {
        min = Q->vec[i];
        pos |= i;
        pos <<= sizeof(unsigned) * 8;
        pos |= P[i];
    }

    /* just to verify correct nj */
    /*
    pos = 0;
    min = DBL_MAX;
    */

    i = D->n;
    Qptr = Q->vec + i;
    while (--i) {
        if (*--Qptr < min) {
            /* check if ith row is min */
            updateQ = DBL_MAX;
            mj = 0;
            N_i = N[i];
            Nj = N - 1;
            sD_i = sD->vec[i];
            sDj = sD->vec - 1;
            Dptr = 0;
            Dfptr = 0;
            Dsptr = 0;
            Dbptr = 0;
            if (D->mat) {
                Dptr = D->mat[i] - 1;
            } else if (D->fmat) {
                Dfptr = D->fmat[i] - 1;
            } else if (D->smat) {
                Dsptr = D->smat[i] - 1;
            } else {
                Dbptr = D->bmat[i] - 1;
            }
            j = -1;
            while (++j < i) {
                q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
                if (0 <= q) {
                    q = ((N_i + *++Nj - 4) >> 1) * q - sD_i - *++sDj;
                    if (q <= updateQ) {
                        updateQ = q;
                        mj = j;
                    }
                } else {
                    ++Nj;
                    ++sDj;
                }
            }

            /* save new Q for ith row */
            P[i] = mj;
            *Qptr = updateQ;
            if (updateQ < min) {
                min = updateQ;
                pos = 0;
                pos |= i;
                pos <<= sizeof(unsigned) * 8;
                pos |= mj;
            }
        }
    }

    return pos;
}

long unsigned maxQpair(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i) {

    int j, mj, N_i, *Nj;
    long unsigned pos;
    double max, q, updateQ, sD_i, *sDj, *Dptr, *Qptr;
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
    pos = 0;
    max = -DBL_MAX;
    if (i && max != Q->vec[i]) {
        max = Q->vec[i];
        pos |= i;
        pos <<= sizeof(unsigned) * 8;
        pos |= P[i];
    }

    /* just to verify correct nj */
    /*
    pos = 0;
    max = -DBL_MAX;
    */

    i = D->n;
    Qptr = Q->vec + i;
    while (--i) {
        if (max < *--Qptr) {
            /* check if ith row is min */
            updateQ = -DBL_MAX;
            mj = 0;
            N_i = N[i];
            Nj = N - 1;
            sD_i = sD->vec[i];
            sDj = sD->vec - 1;
            Dptr = 0;
            Dfptr = 0;
            Dsptr = 0;
            Dbptr = 0;
            if (D->mat) {
                Dptr = D->mat[i] - 1;
            } else if (D->fmat) {
                Dfptr = D->fmat[i] - 1;
            } else if (D->smat) {
                Dsptr = D->smat[i] - 1;
            } else {
                Dbptr = D->bmat[i] - 1;
            }
            j = -1;
            while (++j < i) {
                q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
                if (0 <= q) {
                    q = ((N_i + *++Nj - 4) >> 1) * q - sD_i - *++sDj;
                    if (updateQ <= q) {
                        updateQ = q;
                        mj = j;
                    }
                } else {
                    ++Nj;
                    ++sDj;
                }
            }

            /* save new Q for ith row */
            P[i] = mj;
            *Qptr = updateQ;
            if (max < updateQ) {
                max = updateQ;
                pos = 0;
                pos |= i;
                pos <<= sizeof(unsigned) * 8;
                pos |= mj;
            }
        }
    }

    return pos;
}

long unsigned UPGMApair(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i) {

    int j, mj;
    long unsigned pos;
    double min, q, updateQ, *Dptr, *Qptr;
    float *Dfptr;
    short unsigned *Dsptr;
    unsigned char *Dbptr;

    /* init */
    pos = 0;
    min = DBL_MAX;
    if (i && min != Q->vec[i]) {
        min = Q->vec[i];
        pos |= i;
        pos <<= sizeof(unsigned) * 8;
        pos |= P[i];
    }

    /* just to verify correct upgma */
    /*
    pos = 0;
    min = DBL_MAX;
    */

    i = D->n;
    Qptr = Q->vec + i;
    while (--i) {
        if (*--Qptr < min) {
            if (P[i] < 0) {
                /* check if ith row is min */
                updateQ = DBL_MAX;
                mj = 0;
                Dptr = 0;
                Dfptr = 0;
                Dsptr = 0;
                Dbptr = 0;
                if (D->mat) {
                    Dptr = D->mat[i] - 1;
                } else if (D->fmat) {
                    Dfptr = D->fmat[i] - 1;
                } else if (D->smat) {
                    Dsptr = D->smat[i] - 1;
                } else {
                    Dbptr = D->bmat[i] - 1;
                }
                j = -1;
                while (++j < i) {
                    q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
                    if (0 <= q && q <= updateQ) {
                        updateQ = q;
                        mj = j;
                    }
                }

                /* save new Q for ith row */
                P[i] = mj;
                *Qptr = updateQ;
                if (updateQ < min) {
                    min = updateQ;
                    pos = 0;
                    pos |= i;
                    pos <<= sizeof(unsigned) * 8;
                    pos |= mj;
                }
            } else {
                min = *Qptr;
                pos = 0;
                pos |= i;
                pos <<= sizeof(unsigned) * 8;
                pos |= P[i];
            }
        }
    }

    return pos;
}

int nextQminRow(int i, double min, double *Q) {

    /* get next row */
    i = i ? i : 1;
    Q += i;
    while (--i && min < *--Q);

    return i;
}

int nextQmaxRow(int i, double max, double *Q) {

    /* get next row */
    i = i ? i : 1;
    Q += i;
    while (--i && *--Q < max);

    return i;
}

long unsigned minQrow(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, double *min, int i) {

    int j, mj, N_i;
    long unsigned pos;
    double q, Q_i, sD_i, *sDj, *Dptr;
    float *Dfptr;
    short unsigned *Dsptr;
    unsigned char *Dbptr;

    /* init */
    pos = 0;

    /* check if ith row is min */
    if (Q->vec[i] < *min) {
        N_i = N[i];
        --N;
        sD_i = sD->vec[i];
        sDj = sD->vec - 1;
        Dptr = 0;
        Dfptr = 0;
        Dsptr = 0;
        Dbptr = 0;
        if (D->mat) {
            Dptr = D->mat[i] - 1;
        } else if (D->fmat) {
            Dfptr = D->fmat[i] - 1;
        } else if (D->smat) {
            Dsptr = D->smat[i] - 1;
        } else {
            Dbptr = D->bmat[i] - 1;
        }
        Q_i = DBL_MAX;
        mj = 0;
        j = -1;
        while (++j < i) {
            q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
            if (0 <= q) {
                q = ((N_i + *++N - 4) >> 1) * q - sD_i - *++sDj;
                if (q <= Q_i) {
                    Q_i = q;
                    mj = j;
                }
            } else {
                ++N;
                ++sDj;
            }
        }

        /* save new Q for ith row */
        P[i] = mj;
        Q->vec[i] = Q_i;
        if (Q_i < *min) {
            *min = Q_i;
            pos |= i;
            pos <<= sizeof(unsigned) * 8;
            pos |= mj;
        }
    }

    return pos;
}

long unsigned maxQrow(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, double *max, int i) {

    int j, mj, N_i;
    long unsigned pos;
    double q, Q_i, sD_i, *sDj, *Dptr;
    float *Dfptr;
    short unsigned *Dsptr;
    unsigned char *Dbptr;

    /* init */
    pos = 0;

    /* check if ith row is max */
    if (*max < Q->vec[i]) {
        N_i = N[i];
        --N;
        sD_i = sD->vec[i];
        sDj = sD->vec - 1;
        Dptr = 0;
        Dfptr = 0;
        Dsptr = 0;
        Dbptr = 0;
        if (D->mat) {
            Dptr = D->mat[i] - 1;
        } else if (D->fmat) {
            Dfptr = D->fmat[i] - 1;
        } else if (D->smat) {
            Dsptr = D->smat[i] - 1;
        } else {
            Dbptr = D->bmat[i] - 1;
        }
        Q_i = -DBL_MAX;
        mj = 0;
        j = -1;
        while (++j < i) {
            q = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
            if (0 <= q) {
                q = ((N_i + *++N - 4) >> 1) * q - sD_i - *++sDj;
                if (Q_i <= q) {
                    Q_i = q;
                    mj = j;
                }
            } else {
                ++N;
                ++sDj;
            }
        }

        /* save new Q for ith row */
        P[i] = mj;
        if (*max < (Q->vec[i] = Q_i)) {
            *max = Q_i;
            pos |= i;
            pos <<= sizeof(unsigned) * 8;
            pos |= mj;
        }
    }

    return pos;
}

long unsigned UPGMArow(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, double *min, int i) {

    int j, mj;
    long unsigned pos;
    double d, Q_i, *Dptr;
    float *Dfptr;
    short unsigned *Dsptr;
    unsigned char *Dbptr;

    /* init */
    pos = 0;

    /* check if ith row is min */
    if (Q->vec[i] < *min) {
        if (P[i] < 0) {
            Dptr = 0;
            Dfptr = 0;
            Dsptr = 0;
            Dbptr = 0;
            if (D->mat) {
                Dptr = D->mat[i] - 1;
            } else if (D->fmat) {
                Dfptr = D->fmat[i] - 1;
            } else if (D->smat) {
                Dsptr = D->smat[i] - 1;
            } else {
                Dbptr = D->bmat[i] - 1;
            }
            Q_i = DBL_MAX;
            mj = -1;
            j = -1;
            while (++j < i) {
                d = Dptr ? *++Dptr : Dfptr ? *++Dfptr : Dsptr ? uctod(*++Dsptr) : uctod(*++Dbptr);
                if (0 <= d && d <= Q_i) {
                    Q_i = d;
                    mj = j;
                }
            }

            /* save new Q for ith row */
            P[i] = mj;
            Q->vec[i] = Q_i;
            if (Q_i < *min) {
                *min = Q_i;
                pos |= i;
                pos <<= sizeof(unsigned) * 8;
                pos |= mj;
            }
        } else {
            *min = Q->vec[i];
            pos |= i;
            pos <<= sizeof(unsigned) * 8;
            pos |= P[i];
        }
    }

    return pos;
}

int minQbool(double min, double Min, long unsigned pos, long unsigned Pos) {
    return (min <= Min && (min < Min || Pos < pos));
}

int maxQbool(double max, double Max, long unsigned pos, long unsigned Pos) {
    return (Max <= max && (Max < max || Pos < pos));
}

void *minQ_thread(void *arg) {

    static volatile int thread_begin = 1, thread_wait, Lock = 0;
    static int thread_num = 1, next_i;
    static long unsigned Pos;
    static double Min;
    volatile int *lock = &Lock;
    NJthread *thread = arg;
    int i, *N, *P;
    long unsigned pos;
    double min;
    Matrix *D;
    Vector *sD, *Q;

    /*
    Q(i,j) = (n(i) + n(j) - 4) / 2 * D(i,j) - summa(d(i,k)) - summa(d(k,j));
    sD(i) = summa(d(i,k));
    org:
    Q(i,j) = (n - 2) * D(i,j) - summa(d(i,k)) - summa(d(k,j));
    */

    /* init */
    if (!thread) {
        wait_atomic((thread_begin != thread_num));
        thread_num = 0;
        thread_begin = 0;
        return NULL;
    }
    D = thread->D;
    sD = thread->sD;
    Q = thread->Q;
    N = thread->N;
    P = thread->P;

    if (thread->min) {
        lockTime(lock, 2);
        wait_atomic((thread_begin != thread_num));
        /* init start */
        pos = 0;
        min = Qbool == &minQbool ? DBL_MAX : -DBL_MAX;
        if (thread->mj && min != Q->vec[thread->mj]) {
            min = Q->vec[thread->mj];
            pos |= thread->mj;
            pos <<= sizeof(unsigned) * 8;
            pos |= P[thread->mj];
        }
        Pos = pos;
        Min = min;
        next_i = D->n;
        thread_num = thread->num;
        thread_wait = thread_num;
        thread_begin = 0;
        unlock(lock);
    }

    do {
        wait_atomic(thread_begin);
        /* terminate */
        if (!thread_num) {
            return NULL;
        }

        /* init */
        i = D->n;
        pos = 0;
        min = Qbool == &minQbool ? DBL_MAX : -DBL_MAX;

        /* fill in chunks */
        while (i) {
            /* get next chunk */
            lockTime(lock, 2);
            /* check last min */
            if (pos && Qbool(min, Min, pos, Pos)) {
                Min = min;
                Pos = pos;
            } else {
                min = Min;
            }
            /* get next row */
            i = nextQrow(next_i, min, Q->vec);
            next_i = i;
            unlock(lock);

            /* get row Q */
            if (i) {
                pos = Qrow(D, sD, Q, N, P, &min, i);
            }
        }

        /* wait for remaining threads */
        __sync_sub_and_fetch(&thread_wait, 1);
        wait_atomic(thread_wait);
        __sync_add_and_fetch(&thread_begin, 1);
    } while (!thread->min);

    thread->min = Min;
    thread->pos = Pos;

    return NULL;
}

int updateDNJ(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj) {

    int k, p, Dn, N_j;
    double q, min, sD_j, *D_j, **Dmat, *sDvec, *Qvec;
    float *Df_j, **Dfmat;
    short unsigned *Ds_j, **Dsmat;
    unsigned char *Db_j, **Dbmat;

    /* update D, sD and N */
    updateD(D, sD, N, i, j, Li, Lj);

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
    if (D->mat) {
        Dmat = D->mat;
        D_j = Dmat[j] - 1;
    } else if (D->fmat) {
        Dfmat = D->fmat;
        Df_j = Dfmat[j] - 1;
    } else if (D->smat) {
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
    for (k = 0; k < j; ++k) {
        q = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
        if (0 <= q) {
            q = ((N_j + *++N - 4) >> 1) * q - sD_j - *++sDvec;
            if (q <= *Qvec) {
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
    while (Dn != D->n) {
        if (k == Dn) {
            /* skip ith row */
            Dn = D->n;
            ++Qvec;
            ++P;
            ++sDvec;
            ++N;
        }
        while (++k < Dn) {
            q = Dmat ? Dmat[k][j] : Dfmat ? Dfmat[k][j] : Dsmat ? uctod(Dsmat[k][j]) : uctod(Dbmat[k][j]);
            if (0 <= q) {
                q = ((N_j + *++N - 4) >> 1) * q - sD_j - *++sDvec;
                /* check if new node is better */
                if (q <= *++Qvec) {
                    *Qvec = q;
                    *++P = j;
                    if (q <= min) {
                        min = q;
                        p = k;
                    }
                    /*
                    if(*++P < j || q < *Qvec) {
                        *Qvec = q;
                        *P = j;
                    }
                    */
                } else {
                    ++P;
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

int updateDMN(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int i, int j, double Li, double Lj) {

    int k, p, Dn, N_j;
    double q, max, sD_j, *D_j, **Dmat, *sDvec, *Qvec;
    float *Df_j, **Dfmat;
    short unsigned *Ds_j, **Dsmat;
    unsigned char *Db_j, **Dbmat;

    /* update D, sD and N */
    updateD(D, sD, N, i, j, Li, Lj);

    /* update on (first) row */
    Qvec = Q->vec + j;
    *Qvec = -DBL_MAX;
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
    if (D->mat) {
        Dmat = D->mat;
        D_j = Dmat[j] - 1;
    } else if (D->fmat) {
        Dfmat = D->fmat;
        Df_j = Dfmat[j] - 1;
    } else if (D->smat) {
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
    for (k = 0; k < j; ++k) {
        q = D_j ? *++D_j : Df_j ? *++Df_j : Ds_j ? uctod(*++Ds_j) : uctod(*++Db_j);
        if (0 <= q) {
            q = ((N_j + *++N - 4) >> 1) * q - sD_j - *++sDvec;
            if (*Qvec <= q) {
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
    while (Dn != D->n) {
        if (k == Dn) {
            /* skip ith row */
            Dn = D->n;
            ++Qvec;
            ++P;
            ++sDvec;
            ++N;
        }
        while (++k < Dn) {
            q = Dmat ? Dmat[k][j] : Dfmat ? Dfmat[k][j] : Dsmat ? uctod(Dsmat[k][j]) : uctod(Dbmat[k][j]);
            if (0 <= q) {
                q = ((N_j + *++N - 4) >> 1) * q - sD_j - *++sDvec;
                /* check if new node is better */
                if (*++Qvec <= q) {
                    *Qvec = q;
                    *P = j;
                    if (max <= q) {
                        q = max;
                        p = k;
                    }
                    /*
                    if(*++P < j || *Qvec < q) {
                        *Qvec = q;
                        *P = j;
                    }
                    */
                } else {
                    ++P;
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

int DNJ_popArrange(Matrix *D, Vector *sD, Vector *Q, int *N, int *P, int pos) {

    int i, n, p, N_i, *Nptr;
    double q, min, sD_i, *dest, *src, **mat, *sDvec, *Qptr;
    float *fdest, *fsrc, **fmat;
    short unsigned *sdest, *ssrc, **smat;
    unsigned char *bdest, *bsrc, **bmat;

    n = --D->n;
    --sD->n;
    --Q->n;
    if (pos != n) {
        sD->vec[pos] = sD->vec[n];
        N[pos] = N[n];
        Qptr = Q->vec + pos;
        *Qptr = DBL_MAX;
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
        if (D->mat) {
            mat = D->mat;
            /* row to be emptied */
            dest = mat[pos] - 1;
            /* row to be moved up */
            src = mat[n] - 1;
        } else if (D->fmat) {
            fmat = D->fmat;
            /* row to be emptied */
            fdest = fmat[pos] - 1;
            /* row to be moved up */
            fsrc = fmat[n] - 1;
        } else if (D->smat) {
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
        while (++i < pos) {
            q = dest ? (*++dest = *++src) : fdest ? (*++fdest = *++fsrc) : sdest ? uctod((*++sdest = *++ssrc)) : uctod(
                    (*++bdest = *++bsrc));
            if (0 <= q) {
                q = q * ((N_i + *++Nptr - 4) >> 1) - sD_i - *++sDvec;
                if (q <= *Qptr) {
                    *Qptr = q;
                    *P = i;
                }
            } else {
                ++sDvec;
                ++Nptr;
            }
        }

        /* save min */
        p = pos;
        min = *Qptr;

        /* skip pos, as that is removied */
        if (src) {
            ++src;
        } else if (fsrc) {
            ++fsrc;
        } else if (ssrc) {
            ++ssrc;
        } else {
            ++bsrc;
        }
        ++sDvec;
        ++Nptr;

        /* tilt remaining part of last row into column "pos" */
        while (++i < n) {
            q = mat ? (mat[i][pos] = *++src) : fmat ? (fmat[i][pos] = *++fsrc) : smat ? uctod((smat[i][pos] = *++ssrc))
                                                                                      : uctod((bmat[i][pos] = *++bsrc));
            if (0 <= q) {
                q = q * ((N_i + *++Nptr - 4) >> 1) - sD_i - *++sDvec;
                if (q <= *++Qptr) {
                    *Qptr = q;
                    *++P = pos;
                    if (q <= min) {
                        min = q;
                        p = i;
                    }
                    /*
                    if(*++P < pos || q < *Qptr) {
                        *Qptr = q;
                        *P = pos;
                    }
                    */
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

    /* updates Q */
    /* does not seem to do anything.

    remember to initialize Ptr above.
    double *Ptr = P;
    */
    /*
    mat = D->mat;
    fmat = D->fmat;
    smat = D->smat;
    bmat = D->bmat;
    Qptr = Q->vec;
    sDvec = sD->vec;
    min = *sDvec;
    Nptr = N;
    N_i = *Nptr;
    Ptr = P;
    i = D->n;
    while(--i) {
        q = mat ? (*++mat)[*++Ptr] : fmat ? (*++fmat)[*++Ptr] : smat ? uctod((*++smat)[*++Ptr]) : uctod((*++bmat)[*++Ptr]);
        q = q * ((N_i + *++Nptr - 4) >> 1) - *++sDvec - min;
        if(*++Qptr < q) {
            *Qptr = q;
        }
        if(min < *sDvec) {
            min = *sDvec;
            N_i = *Nptr;
        }
    }
    */

    return p;
}

int minPos(double *Q, int i, int j) {
    return (Q[j] < Q[i] || (i < j && Q[j] == Q[i])) ? j : i;
}

int maxPos(double *Q, int i, int j) {
    return (Q[i] < Q[j] || (i < j && Q[i] == Q[j])) ? j : i;
}

int *dnj(Matrix *D, Vector *sD, Vector *Q, int *N, Qseqs **names) {

    int i, j, mi, mj, mask, shift, *P;
    long unsigned pos;
    double Li, Lj;
    Qseqs *tmp;

    /* init */
    N = initDsDQN(D, sD, Q, N);
    P = N + D->n;
    mask = UINT_MAX;
    shift = 8 * sizeof(unsigned);
    pos = pairQ(Q, P);
    j = pos >> shift;

    /*
    minQpair
    initHNJ
    minQ
    updateDNJ
    DNJ_popArrange
    minPos
    */
    /* get pairs */
    while (D->n != 2 && (pos = Qpair(D, sD, Q, N, P, j))) {
        j = pos & mask;
        i = pos >> shift;

        /* get limbs */
        limbLengthPtr(&Li, &Lj, i, j, sD, N,
                      (D->mat ? D->mat[i][j] : D->fmat ? D->fmat[i][j] : D->smat ? uctod(D->smat[i][j]) : uctod(
                              D->bmat[i][j])));

        /* form leaf */
        formNode(names[j], names[i], Lj, Li);

        /* update D and vectors */
        mi = updateDsDQNPtr(D, sD, Q, N, P, i, j, Li, Lj);

        /* rearrange */
        mj = popArrangePtr(D, sD, Q, N, P, i);
        exchange(names[i], names[D->n], tmp);

        if (mj == D->n) {
            j = mi;
        } else if (mi == D->n) {
            j = mj;
        } else {
            j = qPos(Q->vec, mi, mj);
        }
        ltdMatrixShrink(D, D->size - 1);
    }

    if (D->n == 2) {
        formLastNodePtr(*names, names[1],
                        (D->mat ? **(D->mat) : D->fmat ? **(D->fmat) : D->smat ? uctod(**(D->smat)) : uctod(
                                **(D->bmat))));
    } else {
        /* form remaining nodes with undefined distance */
        while (D->n != 1) {
            /* form leaf */
            formLastNodePtr(*names, names[--D->n], -1.0);
        }
    }

    /* verify valid newick format */
    if (*((*names)->seq) != '(') {
        byteshift((*names)->seq, (*names)->len, '(');
    }

    return N;
}

int *dnj_thread(Matrix *D, Vector *sD, Vector *Q, int *N, Qseqs **names, int thread_num) {
    int i, j, shift, *P;
    long unsigned mask, pos;
    double Li, Lj;
    NJthread *threads, *thread;
    Qseqs *tmp;
    if (thread_num == 1) {
        return dnj(D, sD, Q, N, names);
    }
    /* init */
    N = initDsDQN(D, sD, Q, N);
    P = N + D->n;
    mask = UINT_MAX;
    shift = 8 * sizeof(unsigned);
    pos = pairQ(Q, P);
    i = pos & mask;
    j = pos >> shift;

    /* start threads */
    threads = 0;
    thread = threads;
    i = thread_num;
    while (i) {
        thread = smalloc(sizeof(NJthread));
        thread->D = D;
        thread->sD = sD;
        thread->Q = Q;
        thread->N = N;
        thread->P = P;
        thread->min = 0;
        thread->mi = i;
        thread->mj = j;
        thread->pos = pos;
        thread->num = thread_num;
        thread->next = threads;
        threads = thread;
        if (--i && (errno = pthread_create(&thread->id, NULL, minDist_thread, thread))) {
            fprintf(stderr, "Error: %d (%s)\n", errno, strerror(errno));
            fprintf(stderr, "Will continue with %d threads.\n", thread_num - i);
            i = 0;
        }
    }
    thread->id = 0;
    thread->min = 1;
    /* get pairs */
    while (D->n != 2 && (minDist_thread(thread) || thread->pos)) {
        j = thread->pos & mask;
        i = thread->pos >> shift;
        /* get limbs */
        limbLengthPtr(&Li, &Lj, i, j, sD, N,
                      (D->mat ? D->mat[i][j] : D->fmat ? D->fmat[i][j] : D->smat ? uctod(D->smat[i][j]) : uctod(
                              D->bmat[i][j])));
        /* form leaf */
        formNode(names[j], names[i], Lj, Li);
        /* update D and vectors */
        thread->mi = updateDsDQNPtr(D, sD, Q, N, P, i, j, Li, Lj);
        /* rearrange */
        thread->mj = popArrangePtr(D, sD, Q, N, P, i);
        exchange(names[i], names[D->n], tmp);
        /* get position of min estimate */
        threads->min = 1;
        if (thread->mj == D->n) {
            thread->mj = thread->mi;
        } else if (thread->mi == D->n) {
            thread->mi = thread->mj;
        } else {
            thread->mj = qPos(Q->vec, thread->mi, thread->mj);
        }
        thread->mi = P[thread->mj];
        ltdMatrixShrink(D, D->size - 1);
    }
    if (D->n == 2) {
        formLastNodePtr(*names, names[1],
                        (D->mat ? **(D->mat) : D->fmat ? **(D->fmat) : D->smat ? uctod(**(D->smat)) : uctod(
                                **(D->bmat))));
    } else {
        /* form remaining nodes with undefined distance */
        while (D->n != 1) {
            /* form leaf */
            formLastNodePtr(*names, names[--D->n], -1.0);
        }
    }

    /* verify valid newick format */
    if (*((*names)->seq) != '(') {
        byteshift((*names)->seq, (*names)->len, '(');
    }
    /* collect threads */
    minDist_thread(0);
    thread = threads->next;
    free(threads);
    while (thread) {
        threads = thread->next;
        if ((errno = pthread_join(thread->id, NULL))) {
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
