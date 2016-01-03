#ifndef MAIN_H
#define MAIN_H

#include "encoding.h"

typedef struct {
	unsigned int mismatch_count;
	std::vector<enc_t> subseq_vect;
} result_t;

typedef struct {
	std::vector<enc_t> gene_substrs;
	std::string gene;
} gene_t;

#endif
