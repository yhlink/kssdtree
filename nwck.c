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
#include "dnjheaders/filebuff.h"
#include "dnjheaders/nwck.h"
#include "dnjheaders/qseqs.h"
#include "dnjheaders/pherror.h"

void (*formLastNodePtr)(Qseqs *, Qseqs *, double) = &formLastNode;

void setPrecisionNwck(int precision) {
	
	formNode(0, 0, precision, precision);
	formLastNode(0, 0, precision);
	formLastBiNode(0, 0, precision);
}

void formNode(Qseqs *node1, Qseqs *node2, double L1, double L2) {
	
	static int precision = 9;
	unsigned char *seq, *seq2;
	unsigned newsize;
	double L;
	
	if(!node1 && !node2) { /* set precision */
		precision = L1;
		return;
	} else if(node1->size < node2->size) { /* move largest qseq down */
		exchange(node1->seq, node2->seq, seq);
		exchange(node1->size, node2->size, newsize);
		exchange(node1->len, node2->len, newsize);
		exchange(L1, L2, L);
	}
	
	/* alloc */
	newsize = node1->len + node2->len + 32;
	if(node1->size < newsize) {
		if(!(node1->seq = realloc(node1->seq, newsize))) {
			ERROR();
		}
		node1->size = newsize;
	}
	
	/* shift one byte */
	newsize = ++node1->len;
	seq = node1->seq + newsize;
	seq2 = seq - 1;
	*seq = 0;
	while(--newsize) {
		*--seq = *--seq2;
	}
	*seq2 = '(';
	
	/* form node */
	if(L1 < 0 && L2 < 0) {
		node1->len += sprintf((char *) node1->seq + node1->len, ",%s)", (char *) node2->seq);
	} else {
		node1->len += sprintf((char *) node1->seq + node1->len, ":%.*f,%s:%.*f)", precision, L1, (char *) node2->seq, precision, L2);
	}
}

void formLastNode(Qseqs *node1, Qseqs *node2, double L) {
	
	static int precision = 9;
	unsigned char *seq;
	unsigned newsize;
	
	if(!node1 && !node2) { /* set precision */
		precision = L;
		return;
	} else if(node1->size < node2->size) { /* move largest qseq down */
		exchange(node1->seq, node2->seq, seq);
		exchange(node1->size, node2->size, newsize);
		exchange(node1->len, node2->len, newsize);
	}
	
	/* alloc */
	newsize = node1->len + node2->len + 32;
	if(node1->size < newsize) {
		if(!(node1->seq = realloc(node1->seq, newsize))) {
			ERROR();
		}
		node1->size = newsize;
	}
	
	/* truncate node */
	node1->seq[--node1->len] = 0;
	
	/* form node */
	if(L < 0) {
		node1->len += sprintf((char *) node1->seq + node1->len, ",%s)", (char *) node2->seq);
	} else {
		node1->len += sprintf((char *) node1->seq + node1->len, ",%s:%.*f)", (char *) node2->seq, precision, L);
	}
}

void formLastBiNode(Qseqs *node1, Qseqs *node2, double L) {
	
	static int precision = 9;
	unsigned char *seq, *seq2;
	unsigned newsize;
	
	if(!node1 && !node2) { /* set precision */
		precision = L;
		return;
	} else if(node1->size < node2->size) { /* move largest qseq down */
		exchange(node1->seq, node2->seq, seq);
		exchange(node1->size, node2->size, newsize);
		exchange(node1->len, node2->len, newsize);
	}
	
	/* alloc */
	newsize = node1->len + node2->len + 32;
	if(node1->size < newsize) {
		if(!(node1->seq = realloc(node1->seq, newsize))) {
			ERROR();
		}
		node1->size = newsize;
	}
	
	/* shift one byte */
	newsize = ++node1->len;
	seq = node1->seq + newsize;
	seq2 = seq - 1;
	*seq = 0;
	while(--newsize) {
		*--seq = *--seq2;
	}
	*seq2 = '(';
	
	/* form node */
	if(L < 0) {
		node1->len += sprintf((char *) node1->seq + node1->len, ",%s)", (char *) node2->seq);
	} else {
		L /= 2;
		node1->len += sprintf((char *) node1->seq + node1->len, ":%.*f,%s:%.*f)", precision, L, (char *) node2->seq, precision, L);
	}
}

int getNwck(FileBuff *infile, Qseqs *dest, Qseqs *header) {
	
	int avail, size, len, (*buffFileBuff)(FileBuff *);
	unsigned char *buff, *seq;
	
	/* init */
	avail = infile->bytes;
	buff = infile->next;
	buffFileBuff = infile->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(infile)) == 0) {
			return 0;
		}
		buff = infile->buffer;
	}
	
	/* get name of entry */
	size = header->size;
	seq = header->seq - 1;
	while((*++seq = *buff++) != '(') {
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				return 0;
			}
			buff = infile->buffer;
		}
		if(--size == 0) {
			size = header->size;
			header->size <<= 1;
			if(!(header->seq = realloc(header->seq, header->size))) {
				ERROR();
			}
			seq = header->seq + size - 1;
		}
	}
	*seq = 0;
	header->len = header->size - size;
	
	if(--avail == 0) {
		if((avail = buffFileBuff(infile)) == 0) {
			return 0;
		}
		buff = infile->buffer;
	}
	
	/* get entry */
	size = dest->size;
	seq = dest->seq - 1;
	while((*++seq = *buff++) != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				return 0;
			}
			buff = infile->buffer;
		}
		if(--size == 0) {
			size = dest->size;
			dest->size <<= 1;
			if(!(dest->seq = realloc(dest->seq, dest->size))) {
				ERROR();
			}
			seq = dest->seq + size - 1;
		}
	}
	len = dest->size - size;
	while(--len && *--seq != ')');
	*seq = 0;
	dest->len = len;
	
	infile->bytes = avail - 1;
	infile->next = buff;
	
	return 1;
}

int getSizeNwck(Qseqs *src) {
	
	int n, len;
	unsigned char *ptr;
	
	n = 1;
	ptr = src->seq - 1;
	len = src->len + 1;
	while(--len) {
		if(*++ptr == ',') {
			++n;
		}
	}
	
	return n;
}

double getLimbNwck(Qseqs *node) {
	
	int len;
	double L;
	char *errorMsg;
	unsigned char *seq;
	
	len = node->len;
	seq = node->seq + len;
	if(len-- == 0 || *--seq == ')') {
		return -1.0;
	}
	
	/* find start of limb length */
	while(--len && *--seq != ':');
	
	
	if(len) {
		/* truncate name */
		*seq = 0;
		node->len = len;
		
		/* get limb length */
		L = strtod((char *) ++seq, &errorMsg);
		if(*errorMsg) {
			fprintf(stderr, "Invalid limb length at node:\t%s\n", node->seq);
			exit(1);
		}
	} else {
		L = -1.0;
	}
	
	return L;
}

int stripNwck(Qseqs *node) {
	
	if(*(node->seq) == '(' && node->seq[node->len - 1] == ')') {
		node->len -= 2;
		++node->seq;
		node->seq[node->len] = 0;
		return node->len;
	}
	
	return 0;
}

int splitNwck(Qseqs *node_i, Qseqs *node_j, double *Li, double *Lj) {
	
	int stop, len;
	unsigned char *seq;
	
	/* init */
	len = node_i->len;
	seq = node_i->seq + len;
	
	/* check if split is possible */
	if(!len) {
		return 0;
	}
	
	/* find start of last sub-node in node */
	stop = 0;
	while(stop <= 0 && 0 <= --len) {
		if(*--seq == ')') {
			--stop;
		} else if(*seq == '(') {
			++stop;
		} else if(*seq == ',' && stop == 0) {
			++stop;
		}
	}
	
	/* possible singleton */
	if(stop == 0) {
		return stripNwck(node_i) && splitNwck(node_i, node_j, Li, Lj);
	}
	
	/* truncate org node and move last sub-node to new node */
	*seq = 0;
	node_j->len = node_i->len - len - 2;
	node_j->seq = seq + 1;
	node_i->len = len;
	
	/* check if node is not bifurcating */
	stop = 0;
	while(stop <= 0 && 0 <= --len) {
		if(*--seq == ')') {
			--stop;
		} else if(*seq == '(') {
			++stop;
		} else if(*seq == ',' && stop == 0) {
			++stop;
		}
	}
	
	if(stop != 0) {
		/* not bifurcating */
		*Li = 0;
		*Lj = getLimbNwck(node_j);
	} else {
		/* get limblengths */
		*Li = getLimbNwck(node_i);
		*Lj = getLimbNwck(node_j);
		if(*Lj < 0 && 0 <= *Li) {
			*Lj = 0;
		}
	}
	
	return 1;
}
