#include <string>
#include <fstream>
#include <vector>

#include "parse_file.h"

std::string read_genome(std::ifstream &samples);

std::vector<std::string> parse_file(const char *filename)
{
	std::vector<std::string> genes;
	std::string gene_s;

	std::ifstream samples("shewanella_phages.fasta");
	while(!((gene_s = read_genome(samples)).empty())) {
		genes.push_back(gene_s);
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
