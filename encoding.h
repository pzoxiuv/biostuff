#ifndef ENCODING_H
#define ENCODING_H

#include <string>

typedef struct {
	uint64_t enc_1, enc_2;
} enc_t;

enc_t enc_substr(std::string s);
std::string deencode(enc_t e);

#endif
