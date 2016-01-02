#ifndef MAIN_H
#define MAIN_H

#include "encoding.h"

typedef struct {
	unsigned int mismatch_count;
	std::vector<enc_t> subseq_vect;
} result_t;

#endif
