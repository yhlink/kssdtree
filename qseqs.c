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

#include <stdlib.h>
#include "dnjheaders/pherror.h"
#include "dnjheaders/qseqs.h"

Qseqs * setQseqs(unsigned size) {
	
	Qseqs *dest;
	
	dest = smalloc(sizeof(Qseqs));
	dest->len = 0;
	dest->size = size;
	dest->seq = smalloc(size);
	
	return dest;
}

void destroyQseqs(Qseqs *dest) {
	
	free(dest->seq);
	free(dest);
}

void insertKmerBound(Qseqs *header, int start, int end) {
	
	int *seq;
	
	if((header->len + 2 * sizeof(int)) < header->size) {
		header->size = (header->len + 2 * sizeof(int)) << 1;
		if(!(header->seq = realloc(header->seq, header->size))) {
			ERROR();
		}
	}
	
	seq = (int *) (header->seq + header->len + 1);
	*seq = start;
	*++seq = end;
	header->len += (2 * sizeof(int));
	
}

void qseq2nibble(Qseqs *src, long unsigned *dest) {
	
	int i, j, end, len;
	long unsigned nuc;
	unsigned char *seq;
	
	len = src->len;
	seq = src->seq - 1;
	--dest;
	for(i = 0; i < len; i += 32) {
		end = (i + 32 < len) ? i + 32 : len;
		nuc = 0;
		for(j = i; j < end; ++j) {
			if(*++seq == 4) {
				nuc <<= 2;
			} else {
				nuc = (nuc << 2) | *seq;
			}
		}
		*++dest = nuc;
	}
	if(len & 31) {
		*dest <<= (64 - ((len & 31) << 1));
	}
}
