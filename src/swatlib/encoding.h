#ifndef ENCODING_H
#define ENCODING_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>
#include "nlohmann/json.hpp"

namespace swatlib {

enum class DataType {
    nucleicAcid,
    aminoAcid,
    string  // general string comparison
};

void to_json(nlohmann::json& j, const DataType& d);

void from_json (const nlohmann::json& j, DataType& d);

DataType strToDataType(std::string);

std::string dataTypeToStr(DataType);

// based on BLAST accepted input formats: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp
const std::map<char, uint8_t> NA_CODE = {
	{'N', 0}, // any
	{'A', 1},
	{'T', 2},
	{'C', 3},
	{'G', 4},
	{'U', 5},
	{'R', 6}, // purine (A/G)
	{'Y', 7}, // pyrimidine (T/C)
	{'K', 8}, // ketone (G/T)
	{'M', 9}, // amino (A/C)
	{'S', 10}, // strong (G/C)
	{'W', 11}, // weak (A/T)
	{'B', 12}, // G/T/C
	{'D', 13}, // G/A/T
	{'H', 14}, // A/C/T
	{'V', 15}, // G/C/A
	{'\0', 16}, // stop
};

const std::map<char, uint8_t> AA_CODE = {
	{'A', 0},
	{'R', 1},
	{'N', 2},
	{'D', 3},
	{'C', 4},
	{'Q', 5},
	{'E', 6},
	{'G', 7},
	{'H', 8},
	{'I', 9},
	{'L', 10},
	{'K', 11},
	{'M', 12},
	{'F', 13},
	{'P', 14},
	{'S', 15},
	{'T', 16},
	{'W', 17},
	{'Y', 18},
	{'V', 19},
	{'B', 20},
	{'Z', 21},
	{'X', 22},
	{'*', 23},
	{'\0', 24}
};

class Encoding {
	std::map<char, uint8_t> code;
	std::vector<uint8_t> codeTable;
	std::vector<char> decodeTable;
public:
	Encoding(std::map<char, uint8_t>);
	uint8_t encode(char);
	std::vector<uint8_t> encode(const std::string&);
	std::vector<std::vector<uint8_t>> encode(const std::vector<std::string>&);
	char decode(uint8_t);
	std::string decode(const std::vector<uint8_t>&);
	std::vector<std::string> decode(const std::vector<std::vector<uint8_t>>&);
	uint8_t get_terminator();

	std::vector<uint8_t> getCodeTable();
};

class NAEncoding : public Encoding {
public:
  NAEncoding();
};
class AAEncoding : public Encoding {
public:
	AAEncoding();
};
class CharEncoding : public Encoding {
public:
	CharEncoding();
};

Encoding getEncoder(DataType);
}

#endif // ENCODING_H