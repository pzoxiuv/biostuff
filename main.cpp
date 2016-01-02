#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <climits>
#include <string.h>

typedef struct {
	uint64_t enc_1, enc_2;
} enc_t;

typedef struct {
	unsigned int mismatch_count;
	//std::vector<std::string> subseq_vect;
	std::vector<enc_t> subseq_vect;
} result_t;

std::string read_genome(std::ifstream &samples);
unsigned int count_mismatches(const std::string& s1, const std::string& s2,
					unsigned int len, unsigned int bail_len);
unsigned int count_mismatches2(const char *s1, const char *s2,
					unsigned int len, unsigned int bail_len);
unsigned int count_mismatches3(enc_t s1, enc_t s2,
				unsigned int len, unsigned int bail_len);
bool results_comp(result_t r1, result_t r2);
void print_results(std::vector<result_t> results_vect);
float calc_entropy(result_t res);
enc_t enc_substr(std::string s);
std::string deencode(enc_t e);

const int MAX_MISMATCHES = 7;

const uint8_t a_enc = 0x1;
const uint8_t c_enc = 0x2;
const uint8_t t_enc = 0x4;
const uint8_t g_enc = 0x8;
const uint8_t a_dec = 'A';
const uint8_t c_dec = 'C';
const uint8_t t_dec = 'T';
const uint8_t g_dec = 'G';

int main(int argc, char **argv)
{
	unsigned int substr_size = 20;
	unsigned int num_genes;
	unsigned int i, j, k;

	//std::string cur_gene;
	std::string gene_s;
	//std::vector<char *> cur_gene;
	std::vector<enc_t >cur_gene;
	std::string genes[39];	// what would be more idiomatic than an array?

	//std::vector<char *> gene_substrs[39];
	std::vector<enc_t> gene_substrs[39];
	std::vector<result_t> results_vect;

	std::ifstream samples("shewanella_phages.fasta");

	num_genes = 0;
	while(!((gene_s = read_genome(samples)).empty())) {
		genes[num_genes] = gene_s;
		num_genes++;
	}

	// iterate through each substr of gene ("sliding window"), add to its list of substrs
	for (i=0; i<num_genes; i++) {
			for (j=0; j<(genes[i].size() - substr_size + 1); j++) {
				//char *s = new char[20];	// need these till the end so won't bother delete'ing them
				//memcpy(s, (genes[i].substr(j, substr_size)).c_str(), substr_size);
				//gene_substrs[i].push_back(s);
				gene_substrs[i].push_back(enc_substr(genes[i].substr(j, substr_size)));
			}
	}

	unsigned int best_match = INT_MAX;

	for (k=0; k<num_genes; k++) {
		cur_gene = gene_substrs[k];

		for (i=0; i<cur_gene.size(); i++) {
			//char *cur_substr = cur_gene[i];

			unsigned int least_mismatches[num_genes];
			//std::string least_mismatches_strs[num_genes];
			enc_t least_mismatches_strs[num_genes];

			// for each gene, find substr in that gene with least mismatches to current subseq in cur_gene
			#pragma omp parallel for
			for (j=0; j<num_genes; j++) {
				if (j == k) continue;

				least_mismatches[j] = INT_MAX;

				// iterate over substrs, see how many mismatches they have with current first gene's substr
				//for (auto it2 = gene_substrs[j].begin(); it2 != gene_substrs[j].end(); ++it2) {
				for (unsigned int m = 0; m < gene_substrs[j].size(); m++) {
					//char *gene_substr = *it2;
					//char *gene_substr = gene_substrs[j][m];
					//unsigned int mismatches = count_mismatches(cur_substr, gene_substr,
					//				substr_size, least_mismatches[j]);
					//unsigned int mismatches = count_mismatches2(cur_substr, gene_substr,
					//unsigned int mismatches = count_mismatches2(cur_gene[i], gene_substrs[j][m],
					//				substr_size, least_mismatches[j]);
					unsigned int mismatches = count_mismatches3(cur_gene[i], gene_substrs[j][m],
									substr_size, least_mismatches[j]);

					if (mismatches < least_mismatches[j]) {
						least_mismatches[j] = mismatches;
						//least_mismatches_strs[j] = gene_substr;
						least_mismatches_strs[j] = gene_substrs[j][m];
					}
				}
			}

			result_t cur_res;
			cur_res.mismatch_count = 0;
			for (j=0; j<num_genes; j++) {
				if (j == k) continue;
				cur_res.mismatch_count += least_mismatches[j];
				cur_res.subseq_vect.push_back(least_mismatches_strs[j]);
			}
			//cur_res.subseq_vect.push_back(cur_substr);
			cur_res.subseq_vect.push_back(cur_gene[i]);

			results_vect.push_back(cur_res);

			// If the current mismatch count is worse than the best mismatch count,
			// we know that the next (best match - current match) subsequences won't
			// be better than our current substring - so skip ahead
			if (cur_res.mismatch_count < best_match) {
				best_match = cur_res.mismatch_count;
			} else {
				i += (cur_res.mismatch_count - best_match);
			}
		}
	}

	print_results(results_vect);

	return 0;
}

std::string read_genome(std::ifstream &samples)
{
	std::string s = "", cur_line;

	while (true) {
		std::getline(samples, cur_line);
		if (cur_line.empty()) {
			break;
		} else if (cur_line.find("genome") != std::string::npos) {
			continue;
		} else {
			s += cur_line;
		}
	}
	return s;
}

unsigned int count_mismatches(const std::string& s1, const std::string& s2,
				unsigned int len, unsigned int bail_len)
{
	unsigned int m = 0, i = 0;
	for (i=0; i<len; i++) {
		if (s1[i] != s2[i]) {
			m++;
			if (m >= bail_len) break;	// change to strictly gt to get all substrs that are best len
		}
	}
	return m;
}

unsigned int count_mismatches2(const char *s1, const char *s2,
				unsigned int len, unsigned int bail_len)
{
	unsigned int m = 0, i = 0;
	for (i=0; i<len; i++) {
		if (s1[i] != s2[i]) {
			m++;
			if (m >= bail_len) break;	// change to strictly gt to get all substrs that are best len
		}
	}
	return m;
}

unsigned int count_mismatches3(enc_t s1, enc_t s2,
				unsigned int len, unsigned int bail_len)
{
	uint64_t cnt1 = 0, cnt2 = 0;
	asm volatile (
			"xor %2, %3	\n\t"
			"xor %4, %5	\n\t"
			"popcnt %3, %0	\n\t"
			"popcnt %5, %1	\n\t"
			: "+r"(cnt1), "+r"(cnt2)
			: "r"(s1.enc_1), "r"(s2.enc_1), "r"(s1.enc_2), "r"(s2.enc_2));

	return (cnt1 + cnt2) >> 1;
}

bool results_comp(result_t r1, result_t r2)
{
	return calc_entropy(r1) < calc_entropy(r2);
}

void print_results(std::vector<result_t> results_vect)
{
	unsigned int i, j;

	std::sort(results_vect.begin(), results_vect.end(), &results_comp);

	for (i=0; i<10; i++) {
		std::cout << calc_entropy(results_vect[i]) << ":\n";
		for (j=0; j<results_vect[i].subseq_vect.size(); j++) {
			//std::cout << "\t" << results_vect[i].subseq_vect[j] << "\n";
			std::cout << "\t" << deencode(results_vect[i].subseq_vect[j]) << "\n";
		}
	}
}

float calc_entropy(result_t res)
{
	/*
	 *	For each position in the subsequences
	 *		For each letter:
	 *			1. Count how many occurances of that letter there are
	 *			2. Divide by the number of subsequences
	 *			3. Multiply that number by the log2 of that number
	 *			4. Add to sum
	 */

	const char alphabet[] = {'A', 'C', 'T', 'G'};
	float ent = 0;
	float tmp;

	std::vector<std::string> deencoded_strs;
	for (unsigned int i=0; i<res.subseq_vect.size(); i++) {
		deencoded_strs.push_back(deencode(res.subseq_vect[i]));
	}

	for (unsigned int i=0; i<deencoded_strs[0].size(); i++) {	// char index into subseq
		for (unsigned int j=0; j<4; j++) {	// current letter
			tmp = 0;
			for (unsigned int k=0; k < res.subseq_vect.size(); k++) {	// current subseq
				//if (res.subseq_vect[k][i] == alphabet[j]) {
				if (deencoded_strs[k][i] == alphabet[j]) {
					tmp++;
				}
			}
			tmp /= res.subseq_vect.size();
			if (tmp > 0) {
				ent += (tmp * log2(tmp));
			}
		}
	}

	return -ent;
}

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
