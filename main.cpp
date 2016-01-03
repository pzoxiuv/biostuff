#include <string>
#include <iostream>
#include <algorithm>
#include <climits>

#include "encoding.h"
#include "entropy.h"
#include "main.h"
#include "parse_file.h"

unsigned int count_mismatches(enc_t s1, enc_t s2);
bool results_comp(result_t r1, result_t r2);
void print_results(std::vector<result_t> results_vect);

int main(int argc, char **argv)
{
	unsigned int i, j, k;
	unsigned int best_match = INT_MAX;

	std::vector<result_t> results_vect;
	std::vector<enc_t >cur_gene;
	std::vector<gene_t> genes = parse_file("shewanella_phages.fasta");

	unsigned int *least_mismatches = new unsigned int[genes.size()];	// needed till end, don't delete
	enc_t *least_mismatches_strs = new enc_t[genes.size()];

	for (i=0; i<genes.size(); i++) {
		cur_gene = genes[i].gene_substrs;

		for (j=0; j<cur_gene.size(); j++) {

			// for each gene, find substr in that gene with least mismatches to current subseq in cur_gene
			#pragma omp parallel for
			for (k=0; k<genes.size(); k++) {
				if (k == i) continue;

				least_mismatches[k] = INT_MAX;

				// iterate over substrs, see how many mismatches they have with current first gene's substr
				for (unsigned int m = 0; m < genes[k].gene_substrs.size(); m++) {
					unsigned int mismatches = count_mismatches(cur_gene[j], genes[k].gene_substrs[m]);

					if (mismatches < least_mismatches[k]) {
						least_mismatches[k] = mismatches;
						least_mismatches_strs[k] = genes[k].gene_substrs[m];
					}
				}
			}

			result_t cur_res;
			cur_res.mismatch_count = 0;
			for (k=0; k<genes.size(); k++) {
				if (k == i) continue;
				cur_res.mismatch_count += least_mismatches[k];
				cur_res.subseq_vect.push_back(least_mismatches_strs[k]);
			}
			cur_res.subseq_vect.push_back(cur_gene[j]);

			results_vect.push_back(cur_res);

			// If the current mismatch count is worse than the best mismatch count,
			// we know that the next (best match - current match) subsequences won't
			// be better than our current substring - so skip ahead
			if (cur_res.mismatch_count < best_match) {
				best_match = cur_res.mismatch_count;
			} else {
				j += (cur_res.mismatch_count - best_match);
			}
		}
	}

	print_results(results_vect);

	return 0;
}

unsigned int count_mismatches(enc_t s1, enc_t s2)
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

bool results_comp(result_t r1, result_t r2)
{
	return calc_entropy(r1) < calc_entropy(r2);
}
