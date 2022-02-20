#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <fstream>
#include <string>

// based on specification in: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
namespace swatlib {

struct Fasta {
    std::string description;
    std::string sequence;
};

Fasta readFasta(const std::string& path);

Fasta readFasta(std::ifstream& rf);

}

std::istream &operator>>(std::istream &input, swatlib::Fasta &f);

#endif // FASTA_H