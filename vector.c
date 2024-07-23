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
#include "dnjheaders/vector.h"

Vector * vector_init(unsigned size) {
	
	Vector *dest;
	
	dest = smalloc(sizeof(Vector));
	dest->n = 0;
	dest->size = size;
	dest->vec = smalloc(size * sizeof(double));
	
	return dest;
}

void vector_realloc(Vector *src, unsigned newsize) {
	
	src->vec = realloc(src->vec, newsize * sizeof(double));
	if(!src->vec) {
		ERROR();
	}
	src->size = newsize;
}

void vector_destroy(Vector *src) {
	
	free(src->vec);
	free(src);
}
