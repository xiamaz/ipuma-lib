#include <algorithm>
#include "fasta.h"

void rTrim(std::string& str) {
    str.erase(
        find_if(str.rbegin(), str.rend(), [](int c) {return !isspace(c);}).base(),
        str.end()
    );
}

void removeWhitespace(std::string& str) {
    str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
}

std::istream& operator>>(std::istream& input, swatlib::Fasta& f)
{
    std::string buffer;
    bool seenHeader = false;
    int lines = 0;
    f.sequence = "";
    int len = input.tellg();
    while (getline(input, buffer)) {
        if (!seenHeader && buffer[0] == '>') {
            seenHeader = true;
            rTrim(buffer);
            f.description = buffer;
        } else if (seenHeader && buffer[0] == '>') {
            input.seekg(len ,std::ios_base::beg);
            break;
        } else if (seenHeader) {
            removeWhitespace(buffer);
            f.sequence += buffer;
        }
        lines++;
        len = input.tellg();
    }
    return input;
}

swatlib::Fasta swatlib::readFasta(std::ifstream& rf) {
    swatlib::Fasta f;
    rf >> f;
    return f;
}

swatlib::Fasta swatlib::readFasta(const std::string& path) {
    std::ifstream rf;
    rf.open(path);
    Fasta f = readFasta(rf);
    rf.close();
    return f;
}