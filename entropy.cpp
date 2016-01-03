#include <algorithm>

#include "entropy.h"
#include "main.h"

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
		for (unsigned int j=0; j<sizeof(alphabet); j++) {	// current letter
			tmp = 0;
			for (unsigned int k=0; k < res.subseq_vect.size(); k++) {	// current subseq
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
