#ifndef IPUMA_TYPES_H
#define IPUMA_TYPES_H

typedef std::vector<std::string> RawSequences;
typedef std::vector<uint8_t> EncSequences;

struct __attribute__((__packed__)) Comparison {
	int32_t indexA;
	int32_t indexB;
};

typedef std::vector<Comparison> Comparisons;

#endif