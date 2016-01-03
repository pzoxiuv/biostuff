#include <string>
#include <iostream>
#include <algorithm>
#include <climits>

#include "encoding.h"
#include "entropy.h"
#include "main.h"
#include "parse_file.h"

unsigned int count_mismatches(enc_t s1, enc_t s2,
				unsigned int len, unsigned int bail_len);
bool results_comp(result_t r1, result_t r2);
void print_results(std::vector<result_t> results_vect);

int main(int argc, char **argv)
{
	unsigned int substr_size = 20;
	unsigned int i, j, k;

	std::vector<enc_t >cur_gene;

	std::vector<enc_t> gene_substrs[39];
	std::vector<result_t> results_vect;

	std::vector<std::string> genes = parse_file("shewanella_phages.fasta");

	// iterate through each substr of gene ("sliding window"), add to its list of substrs
	for (i=0; i<genes.size(); i++) {
			for (j=0; j<(genes[i].size() - substr_size + 1); j++) {
				gene_substrs[i].push_back(enc_substr(genes[i].substr(j, substr_size)));
			}
	}

	unsigned int best_match = INT_MAX;

	for (k=0; k<genes.size(); k++) {
		cur_gene = gene_substrs[k];

		for (i=0; i<cur_gene.size(); i++) {

			unsigned int least_mismatches[genes.size()];
			enc_t least_mismatches_strs[genes.size()];

			// for each gene, find substr in that gene with least mismatches to current subseq in cur_gene
			#pragma omp parallel for
			for (j=0; j<genes.size(); j++) {
				if (j == k) continue;

				least_mismatches[j] = INT_MAX;

				// iterate over substrs, see how many mismatches they have with current first gene's substr
				for (unsigned int m = 0; m < gene_substrs[j].size(); m++) {
					unsigned int mismatches = count_mismatches(cur_gene[i], gene_substrs[j][m],
									substr_size, least_mismatches[j]);

					if (mismatches < least_mismatches[j]) {
						least_mismatches[j] = mismatches;
						least_mismatches_strs[j] = gene_substrs[j][m];
					}
				}
			}

			result_t cur_res;
			cur_res.mismatch_count = 0;
			for (j=0; j<genes.size(); j++) {
				if (j == k) continue;
				cur_res.mismatch_count += least_mismatches[j];
				cur_res.subseq_vect.push_back(least_mismatches_strs[j]);
			}
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

unsigned int count_mismatches(enc_t s1, enc_t s2,
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
			std::cout << "\t" << deencode(results_vect[i].subseq_vect[j]) << "\n";
		}
	}
}
