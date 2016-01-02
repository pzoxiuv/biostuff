#include <string>
#include "encoding.h"

const uint8_t a_enc = 0x1;
const uint8_t c_enc = 0x2;
const uint8_t t_enc = 0x4;
const uint8_t g_enc = 0x8;
const uint8_t a_dec = 'A';
const uint8_t c_dec = 'C';
const uint8_t t_dec = 'T';
const uint8_t g_dec = 'G';

enc_t enc_substr(std::string s)
{
	uint8_t *s_enc = new uint8_t[10];

	for (unsigned int i=0; i<s.size(); i+=2) {
		uint8_t tmp1 = 0, tmp2 = 0;
		if (s[i] == a_dec) tmp1 = a_enc;
		else if (s[i] == c_dec) tmp1 = c_enc;
		else if (s[i] == t_dec) tmp1 = t_enc;
		else if (s[i] == g_dec) tmp1 = g_enc;
		if (s[i+1] == a_dec) tmp2 = a_enc;
		else if (s[i+1] == c_dec) tmp2 = c_enc;
		else if (s[i+1] == t_dec) tmp2 = t_enc;
		else if (s[i+1] == g_dec) tmp2 = g_enc;
		s_enc[i/2] = ((tmp2 & 0x0f) << 4) | (tmp1 & 0x0f);
	}

	uint64_t s_1 = 0;
	for (unsigned int i=0; i<8; i++) {
		s_1 |= ((uint64_t) s_enc[i]) << (i*8);
	}
	uint64_t s_2 = (((uint64_t) s_enc[9]) << 8) | (uint64_t) s_enc[8];

	enc_t encoded;
	encoded.enc_1 = s_1;
	encoded.enc_2 = s_2;

	return encoded;
}

std::string deencode(enc_t e)
{
	std::string res;

	for (unsigned int i=0; i<16; i++) {
		if ((e.enc_1 & 0xf) == a_enc) res += a_dec;
		if ((e.enc_1 & 0xf) == c_enc) res += c_dec;
		if ((e.enc_1 & 0xf) == t_enc) res += t_dec;
		if ((e.enc_1 & 0xf) == g_enc) res += g_dec;
		e.enc_1 >>= 4;
	}
	for (unsigned int i=0; i<4; i++) {
		if ((e.enc_2 & 0xf) == a_enc) res += a_dec;
		if ((e.enc_2 & 0xf) == c_enc) res += c_dec;
		if ((e.enc_2 & 0xf) == t_enc) res += t_dec;
		if ((e.enc_2 & 0xf) == g_enc) res += g_dec;
		e.enc_2 >>= 4;
	}

	return res;
}
