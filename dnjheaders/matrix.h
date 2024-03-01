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

#include <stdio.h>

#ifndef MATRIX
typedef struct matrix Matrix;
struct matrix {
	int n;
	int size;
	double **mat;
	float **fmat;
	short unsigned **smat;
	unsigned char **bmat;
	FILE *file;
};
#define MATRIX 1
#endif

extern Matrix * (*ltdMatrix_init)(int);
extern void (*ltdMatrixShrink)(Matrix *, int);

Matrix * ltdMatrixInit(int size);
Matrix * ltdMatrixMinit(int size);
void ltdMatrix_mrealloc(Matrix *src, int size);
void ltdMatrix_realloc(Matrix *src, int size);
void Matrix_mdestroy(Matrix *src);
void Matrix_destroy(Matrix *src);
void ltdMatrix_popArrange(Matrix *mat, unsigned pos);
int ltdMatrix_add(Matrix *src);
void ltdMatrix_shrink(Matrix *src, int size);
void ltdMatrix_nullInt(Matrix *src, int size);
