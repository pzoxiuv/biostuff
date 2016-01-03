#include <string>
#include <fstream>
#include <vector>

#include "encoding.h"
#include "main.h"
#include "parse_file.h"

std::string read_genome(std::ifstream &samples);

const unsigned int SUBSTR_LEN = 20;

std::vector<gene_t> parse_file(const char *filename)
{
	gene_t *cur_gene;
	std::vector<gene_t> genes;
	std::string gene_s;

	std::ifstream samples("shewanella_phages.fasta");
	while(!((gene_s = read_genome(samples)).empty())) {
		cur_gene = new gene_t;	// this is never delete'd, needed till end of program
		cur_gene->gene = gene_s;

		for (unsigned int i=0; i<(gene_s.size() - SUBSTR_LEN + 1); i++) {
			cur_gene->gene_substrs.push_back(enc_substr(gene_s.substr(i, SUBSTR_LEN)));
		}

		genes.push_back(*cur_gene);
	}

	return genes;
}

std::string read_genome(std::ifstream& samples)
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
