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
#include "filebuff.h"
#include "matrix.h"
#include "qseqs.h"

extern char * (*stripEntry)(char *);

char * stripDir(char *str);
char * noStripDir(char *str);
void setPrecisionPhy(int precision);
void printphy(FILE *outfile, Matrix *src, char **names, unsigned char *include, char *comment, unsigned format);
void printphyUpdate(FILE *outfile, int n, char *name, double *D, unsigned format);
void printfullphy(FILE *outfile, Matrix *src, char **names, unsigned format);
Qseqs ** loadPhy(Matrix *src, Qseqs **names, Qseqs *header, FileBuff *infile, char sep, char quotes);
int getSizePhy(FileBuff *infile);
Qseqs ** getFilenamesPhy(char *path, int n, FileBuff *infile, char sep);
