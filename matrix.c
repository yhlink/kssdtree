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
#include <stdlib.h>
#include <stdio.h>
#include <sys/param.h>
#include "dnjheaders/matrix.h"
#include "dnjheaders/pherror.h"
#include "dnjheaders/tmp.h"
#ifdef _WIN32
#include "kssdheaders/mman.h"
#else
#include <sys/mman.h>
#endif
Matrix * (*ltdMatrix_init)(int) = &ltdMatrixInit;
void (*ltdMatrixShrink)(Matrix *, int) = &ltdMatrix_nullInt;

Matrix * ltdMatrixInit(int size) {

	static long unsigned type = sizeof(double);
	int i;
	long unsigned Size;
	double **ptr, *src;
	float **fptr, *fsrc;
	short unsigned **sptr, *ssrc;
	unsigned char **bptr, *bsrc;
	Matrix *dest;

	if(size < 0) {
		type = -size;
		return 0;
	}

	dest = smalloc(sizeof(Matrix));
	dest->n = 0;
	dest->size = size;
	dest->mat = 0;
	dest->fmat = 0;
	dest->smat = 0;
	dest->bmat = 0;
	Size = size;
	Size *= (size - 1);
	Size >>= 1;
	Size *= type;
	if(type == sizeof(double)) {
		dest->mat = smalloc(dest->size * sizeof(double *));
		*(dest->mat) = smalloc(Size);
	} else if(type == sizeof(float)) {
		dest->fmat = smalloc(dest->size * sizeof(float *));
		*(dest->fmat) = smalloc(Size);
	} else if(type == sizeof(short unsigned)) {
		dest->smat = smalloc(dest->size * sizeof(short unsigned *));
		*(dest->smat) = smalloc(Size);
	} else {
		dest->bmat = smalloc(dest->size * sizeof(unsigned char *));
		*(dest->bmat) = smalloc(Size);
	}
	dest->file = 0;

	/* set matrix rows */
	if(dest->mat) {
		ptr = dest->mat;
		src = *ptr;
		i = 0;
		*ptr = src;
		while(--size) {
			*++ptr = src + i;
			src += i++;
		}
	} else if(dest->fmat) {
		fptr = dest->fmat;
		fsrc = *fptr;
		i = 0;
		*fptr = fsrc;
		while(--size) {
			*++fptr = fsrc + i;
			fsrc += i++;
		}
	} else if(dest->smat) {
		sptr = dest->smat;
		ssrc = *sptr;
		i = 0;
		*sptr = ssrc;
		while(--size) {
			*++sptr = ssrc + i;
			ssrc += i++;
		}
	} else {
		bptr = dest->bmat;
		bsrc = *bptr;
		i = 0;
		*bptr = bsrc;
		while(--size) {
			*++bptr = bsrc + i;
			bsrc += i++;
		}
	}

	return dest;
}

Matrix * ltdMatrixMinit(int size) {

	static long unsigned type = sizeof(double);
	int i;
	long unsigned Size;
	double **ptr, *src;
	float **fptr, *fsrc;
	short unsigned **sptr, *ssrc;
	unsigned char **bptr, *bsrc;
	FILE *tmp;
	Matrix *dest;

	if(size < 0) {
		type = -size;
		return 0;
	}

	dest = smalloc(sizeof(Matrix));
	dest->n = 0;
	dest->size = size;
	dest->mat = 0;
	dest->fmat = 0;
	dest->smat = 0;
	dest->bmat = 0;
	if(type == sizeof(double)) {
		dest->mat = smalloc(size * sizeof(double *));
	} else if(type == sizeof(float)) {
		dest->fmat = smalloc(size * sizeof(float *));
	} else if(type == sizeof(short unsigned)) {
		dest->smat = smalloc(size * sizeof(short unsigned *));
	} else {
		dest->bmat = smalloc(size * sizeof(unsigned char *));
	}
	tmp = tmpF(0);
	dest->file = tmp;
	Size = size;
	Size *= (size - 1);
	Size >>= 1;
	Size *= type;
	if(fseek(tmp, Size - 1, SEEK_SET) || putc(0, tmp) == EOF) {
		ERROR();
	}
	fflush(tmp);
	sfseek(tmp, 0, SEEK_SET);
	if(dest->mat) {
		*(dest->mat) = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(*(dest->mat) == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(*(dest->mat), Size, POSIX_MADV_SEQUENTIAL);

		/* set matrix rows */
		ptr = dest->mat;
		src = *ptr;
		i = 0;
		*ptr = src;
		while(--size) {
			*++ptr = src + i;
			src += i++;
		}
	} else if(dest->fmat) {
		*(dest->fmat) = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(*(dest->fmat) == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(*(dest->fmat), Size, POSIX_MADV_SEQUENTIAL);

		/* set matrix rows */
		fptr = dest->fmat;
		fsrc = *fptr;
		i = 0;
		*fptr = fsrc;
		while(--size) {
			*++fptr = fsrc + i;
			fsrc += i++;
		}
	} else if(dest->smat) {
		*(dest->smat) = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(*(dest->smat) == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(*(dest->smat), Size, POSIX_MADV_SEQUENTIAL);

		/* set matrix rows */
		sptr = dest->smat;
		ssrc = *sptr;
		i = 0;
		*sptr = ssrc;
		while(--size) {
			*++sptr = ssrc + i;
			ssrc += i++;
		}
	} else {
		*(dest->bmat) = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(*(dest->bmat) == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(*(dest->bmat), Size, POSIX_MADV_SEQUENTIAL);

		/* set matrix rows */
		bptr = dest->bmat;
		bsrc = *bptr;
		i = 0;
		*bptr = bsrc;
		while(--size) {
			*++bptr = bsrc + i;
			bsrc += i++;
		}
	}

	return dest;
}

void ltdMatrix_mrealloc(Matrix *src, int size) {

	int i;
	long unsigned type, Size;
	double **ptr, *mat;
	float **fptr, *fmat;
	short unsigned **sptr, *smat;
	unsigned char **bptr, *bmat;
	FILE *tmp;

	/* init */
	type = src->mat ? sizeof(double) : src->fmat ? sizeof(float) : src->smat ? sizeof(short unsigned) : sizeof(unsigned char);

	/* unmap current mapping */
	Size = src->size;
	Size *= (src->size - 1);
	Size >>= 1;
	Size *= type;
	if(src->mat) {
		msync(*(src->mat), Size, MS_SYNC);
		munmap(*(src->mat), Size);
	} else if(src->fmat) {
		msync(*(src->fmat), Size, MS_SYNC);
		munmap(*(src->fmat), Size);
	} else if(src->smat) {
		msync(*(src->smat), Size, MS_SYNC);
		munmap(*(src->smat), Size);
	} else {
		msync(*(src->bmat), Size, MS_SYNC);
		munmap(*(src->bmat), Size);
	}

	/* reallocate file */
	tmp = src->file;
	Size = size;
	Size *= (size - 1);
	Size >>= 1;
	Size *= type;
	if(fseek(tmp, Size - 1, SEEK_SET) || putc(0, tmp) == EOF) {
		ERROR();
	}
	fflush(tmp);
	sfseek(tmp, 0, SEEK_SET);

	/* map new size */
	if(type == sizeof(double)) {
		mat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(mat == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
#ifdef __linux__
        posix_madvise(mat, Size, POSIX_MADV_SEQUENTIAL);
#endif
		ptr = realloc(src->mat, size * sizeof(double *));
		if(!ptr) {
			ERROR();
		}
		src->mat = ptr;
		*(src->mat) = mat;
		src->size = size;

		/* set matrix rows */
		i = 0;
		*ptr = mat;
		while(--size) {
			*++ptr = mat + i;
			mat += i++;
		}
	} else if(type == sizeof(float)) {
		fmat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(fmat == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
#ifdef __linux__
        posix_madvise(fmat, Size, POSIX_MADV_SEQUENTIAL);
#endif
		fptr = realloc(src->fmat, size * sizeof(float *));
		if(!fptr) {
			ERROR();
		}
		src->fmat = fptr;
		*(src->fmat) = fmat;
		src->size = size;

		/* set matrix rows */
		i = 0;
		*fptr = fmat;
		while(--size) {
			*++fptr = fmat + i;
			fmat += i++;
		}
	} else if(type == sizeof(short unsigned)) {
		smat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(smat == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
#ifdef __linux__
        posix_madvise(smat, Size, POSIX_MADV_SEQUENTIAL);
#endif
		sptr = realloc(src->smat, size * sizeof(short unsigned *));
		if(!sptr) {
			ERROR();
		}
		src->smat = sptr;
		*(src->smat) = smat;
		src->size = size;

		/* set matrix rows */
		i = 0;
		*sptr = smat;
		while(--size) {
			*++sptr = smat + i;
			smat += i++;
		}
	} else {
		bmat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(bmat == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
#ifdef __linux__
        posix_madvise(bmat, Size, POSIX_MADV_SEQUENTIAL);
#endif
		bptr = realloc(src->bmat, size * sizeof(unsigned char *));
		if(!bptr) {
			ERROR();
		}
		src->bmat = bptr;
		*(src->bmat) = bmat;
		src->size = size;
		
		/* set matrix rows */
		i = 0;
		*bptr = bmat;
		while(--size) {
			*++bptr = bmat + i;
			bmat += i++;
		}
	}
}

void ltdMatrix_realloc(Matrix *src, int size) {
	
	int i;
	long unsigned Size;
	double **ptr, *mat;
	float **fptr, *fmat;
	short unsigned **sptr, *smat;
	unsigned char **bptr, *bmat;
	
	if(src->file) {
		ltdMatrix_mrealloc(src, size);
		return;
	}
	
	Size = size;
	Size *= (size - 1);
	Size >>= 1;
	if(src->mat) {
		Size *= sizeof(double);
		mat = realloc(*(src->mat), Size);
		ptr = realloc(src->mat, size * sizeof(double *));
		if(!ptr || !mat) {
			ERROR();
		}
		src->mat = ptr;
		*(src->mat) = mat;
		src->size = size;
		
		/* set matrix rows */
		i = 0;
		*ptr = mat;
		while(--size) {
			*++ptr = mat + i;
			mat += i++;
		}
	} else if(src->fmat) {
		Size *= sizeof(float);
		fmat = realloc(*(src->fmat), Size);
		fptr = realloc(src->fmat, size * sizeof(float *));
		if(!fptr || !fmat) {
			ERROR();
		}
		src->fmat = fptr;
		*(src->fmat) = fmat;
		src->size = size;
		
		/* set matrix rows */
		i = 0;
		*fptr = fmat;
		while(--size) {
			*++fptr = fmat + i;
			fmat += i++;
		}
	} else if(src->smat) {
		Size *= sizeof(short unsigned);
		smat = realloc(*(src->smat), Size);
		sptr = realloc(src->smat, size * sizeof(short unsigned *));
		if(!sptr || !smat) {
			ERROR();
		}
		src->smat = sptr;
		*(src->smat) = smat;
		src->size = size;
		
		/* set matrix rows */
		i = 0;
		*sptr = smat;
		while(--size) {
			*++sptr = smat + i;
			smat += i++;
		}
	} else {
		Size *= sizeof(unsigned char);
		bmat = realloc(*(src->bmat), Size);
		bptr = realloc(src->bmat, size * sizeof(unsigned char *));
		if(!bptr || !bmat) {
			ERROR();
		}
		src->bmat = bptr;
		*(src->bmat) = bmat;
		src->size = size;
		
		/* set matrix rows */
		i = 0;
		*bptr = bmat;
		while(--size) {
			*++bptr = bmat + i;
			bmat += i++;
		}
	}
}

void Matrix_mdestroy(Matrix *src) {
	
	long unsigned size;
	
	if(src) {
		size = src->size;
		size *= (src->size - 1);
		size >>= 1;
		if(src->mat) {
			size *= sizeof(double);
			msync(*(src->mat), size, MS_SYNC);
			munmap(*(src->mat), size);
			free(src->mat);
		} else if(src->fmat) {
			size *= sizeof(float);
			msync(*(src->fmat), size, MS_SYNC);
			munmap(*(src->fmat), size);
			free(src->fmat);
		} else if(src->smat) {
			size *= sizeof(short unsigned);
			msync(*(src->smat), size, MS_SYNC);
			munmap(*(src->smat), size);
			free(src->smat);
		} else if(src->bmat) {
			size *= sizeof(unsigned char);
			msync(*(src->bmat), size, MS_SYNC);
			munmap(*(src->bmat), size);
			free(src->bmat);
		}
		fclose(src->file);
		free(src);
	}
}

void Matrix_destroy(Matrix *src) {
	
	if(src) {
		if(src->file) {
			Matrix_mdestroy(src);
			return;
		} else if(src->mat) {
			free(*(src->mat));
			free(src->mat);
		} else if(src->fmat) {
			free(*(src->fmat));
			free(src->fmat);
		} else if(src->smat) {
			free(*(src->smat));
			free(src->smat);
		} else if(src->bmat) {
			free(*(src->bmat));
			free(src->bmat);
		}
		free(src);
	}
}

void ltdMatrix_popArrange(Matrix *mat, unsigned pos) {
	
	int i, n;
	double *dest, *src;
	float *fdest, *fsrc;
	short unsigned *sdest, *ssrc;
	unsigned char *bdest, *bsrc;
	
	n = --mat->n;
	if(pos != n) {
		if(mat->mat) {
			/* row to be emptied */
			dest = mat->mat[pos] - 1;
			/* row to be moved up */
			src = mat->mat[n] - 1;
			
			/* copy last row into "first" row */
			i = pos + 1;
			while(--i) {
				*++dest = *++src;
			}
			
			/* tilt remaining part of last row into column "pos" */
			i = pos;
			++src;
			while(++i < n) {
				mat->mat[i][pos] = *++src;
			}
		} else if(mat->fmat) {
			/* row to be emptied */
			fdest = mat->fmat[pos] - 1;
			/* row to be moved up */
			fsrc = mat->fmat[n] - 1;
			
			/* copy last row into "first" row */
			i = pos + 1;
			while(--i) {
				*++fdest = *++fsrc;
			}
			
			/* tilt remaining part of last row into column "pos" */
			i = pos;
			++fsrc;
			while(++i < n) {
				mat->fmat[i][pos] = *++fsrc;
			}
		} else if(mat->smat) {
			/* row to be emptied */
			sdest = mat->smat[pos] - 1;
			/* row to be moved up */
			ssrc = mat->smat[n] - 1;
			
			/* copy last row into "first" row */
			i = pos + 1;
			while(--i) {
				*++sdest = *++ssrc;
			}
			
			/* tilt remaining part of last row into column "pos" */
			i = pos;
			++ssrc;
			while(++i < n) {
				mat->smat[i][pos] = *++ssrc;
			}
		} else {
			/* row to be emptied */
			bdest = mat->bmat[pos] - 1;
			/* row to be moved up */
			bsrc = mat->bmat[n] - 1;
			
			/* copy last row into "first" row */
			i = pos + 1;
			while(--i) {
				*++bdest = *++bsrc;
			}
			
			/* tilt remaining part of last row into column "pos" */
			i = pos;
			++bsrc;
			while(++i < n) {
				mat->bmat[i][pos] = *++bsrc;
			}
		}
	}
}

int ltdMatrix_add(Matrix *src) {
	
	int i;
	double *ptr;
	float *fptr;
	short unsigned *sptr;
	unsigned char *bptr;
	
	/* realloc */
	if(++src->n == src->size) {
		ltdMatrix_realloc(src, src->size << 1);
	}
	
	/* init new row */
	i = src->n;
	if(src->mat) {
		ptr = src->mat[i - 1] - 1;
		while(--i) {
			*++ptr = 0;
		}
	} else if(src->fmat) {
		fptr = src->fmat[i - 1] - 1;
		while(--i) {
			*++fptr = 0;
		}
	} else if(src->smat) {
		sptr = src->smat[i - 1] - 1;
		while(--i) {
			*++sptr = 0;
		}
	} else {
		bptr = src->bmat[i - 1] - 1;
		while(--i) {
			*++bptr = 0;
		}
	}
	
	return src->n - 1;
}

void ltdMatrix_shrink(Matrix *src, int size) {
	
	int i;
	long unsigned Size, type;
	double **ptr, *mat;
	float **fptr, *fmat;
	short unsigned **sptr, *smat;
	unsigned char **bptr, *bmat;
	
	if((size & 2047) || src->size < size) {
		return;
	} else if(!src->file) {
		ltdMatrix_realloc(src, size);
		return;
	}
	
	/* unmap mmap */
	Size = src->size;
	Size *= (src->size - 1);
	Size >>= 1;
	type = src->mat ? sizeof(double) : src->fmat ? sizeof(float) : src->smat ? sizeof(short unsigned) : sizeof(unsigned char);
	Size *= type;
	if(src->mat) {
		ptr = src->mat;
	} else if(src->fmat) {
		ptr = (double **)(src->fmat);
	} else if(src->smat) {
		ptr = (double **)(src->smat);
	} else {
		ptr = (double **)(src->bmat);
	}
	mat = *ptr;
	msync(mat, Size, MS_SYNC);
	munmap(mat, Size);
	
	/* map new size */
	Size = size;
	Size *= (size - 1);
	Size >>= 1;
	Size *= type;
	if(ftruncate(fileno(src->file), Size)) {
		ERROR();
	}
	mat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(src->file), 0);
	if(mat == MAP_FAILED) {
		fprintf(stderr, "MMAP failed:\n");
		ERROR();
	}
	//posix_madvise(mat, Size, POSIX_MADV_SEQUENTIAL);
	if(!(ptr = realloc(ptr, size * sizeof(double *)))) {
		ERROR();
	}
	src->size = size;
	
	/* set matrix rows */
	if(src->mat) {
		src->mat = ptr;
		*(src->mat) = mat;
		i = 0;
		*ptr = mat;
		while(--size) {
			*++ptr = mat + i;
			mat += i++;
		}
	} else if(src->fmat) {
		fptr = (float **)(ptr);
		fmat = (float *)(mat);
		src->fmat = fptr;
		*(src->fmat) = fmat;
		i = 0;
		*fptr = fmat;
		while(--size) {
			*++fptr = fmat + i;
			fmat += i++;
		}
	} else if(src->smat) {
		sptr = (short unsigned **)(ptr);
		smat = (short unsigned *)(mat);
		src->smat = sptr;
		*(src->smat) = smat;
		i = 0;
		*sptr = smat;
		while(--size) {
			*++sptr = smat + i;
			smat += i++;
		}
	} else {
		bptr = (unsigned char **)(ptr);
		bmat = (unsigned char *)(mat);
		src->bmat = bptr;
		*(src->bmat) = bmat;
		i = 0;
		*bptr = bmat;
		while(--size) {
			*++bptr = bmat + i;
			bmat += i++;
		}
	}
}

void ltdMatrix_nullInt(Matrix *src, int size) {}
