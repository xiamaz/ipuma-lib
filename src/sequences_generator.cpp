#include "sequences_generator.h"

#include <random>

namespace ipu {

std::map<swatlib::DataType, std::string> SYMBOLS = {
	{swatlib::DataType::nucleicAcid, "ATCG"},
	{swatlib::DataType::aminoAcid, "ARNDCQEGHILKMFPSTWYVBZX"},
};

std::string generateRandomSequence(size_t seqLen, swatlib::DataType dtype, std::mt19937& gen) {
	const auto& chars = SYMBOLS[dtype];
	std::uniform_int_distribution<> distrib(0, chars.size() - 1);

	std::stringstream ss;
	for (int i = 0; i < seqLen; ++i) {
		ss << chars[distrib(gen)];
	}
	return ss.str();
}

char get_other_symbol(char c, const std::string& symbols) {
	auto i = symbols.find(c);
	if (i == std::string::npos) {
		throw std::runtime_error("Illegal character: (" + c + std::string(":") + std::to_string(c) + std::string(")"));
	}
	return symbols[(i + 1) % symbols.size()];
}

std::string mutateSequence(std::string seq, swatlib::DataType dtype, std::mt19937& gen, float similarity, int seedBegin, int seedSize) {
	std::uniform_real_distribution<float> prob(0, 1);
	const auto& chars = SYMBOLS[dtype];
	for (int i = 0; i < seq.size(); ++i) {
		if (i >= seedBegin && i <= seedBegin + seedSize) continue;
		if (prob(gen) > similarity) {
			seq[i] = get_other_symbol(seq[i], chars);
		}
	}
	return seq;
}

template<typename C>
SequenceDatabase<C> loadSequences(GeneratorConfig& config, SWConfig& swconfig) {
  std::mt19937 gen(config.generatorSeed);

	SequenceDatabase<C> db;

	int seed_start = config.generatorSeqLen / 2 - swconfig.seedLength / 2;
	if (seed_start + swconfig.seedLength >= config.generatorSeqLen) {
		throw std::runtime_error("Seed does not fit into sequence length");
	}
	for (int i = 0; i < config.generatorCount; ++i) {
		std::string seq1, seq2;
		if (config.generatorSimilarity > 0) {
			seq1 = generateRandomSequence(config.generatorSeqLen, swconfig.datatype, gen);
			seq2 = mutateSequence(seq1, swconfig.datatype, gen, config.generatorSimilarity, seed_start, swconfig.seedLength);
		} else {
			const auto& chars = SYMBOLS[swconfig.datatype];
			seq1 = std::string(config.generatorSeqLen, chars[0]);
			seq2 = std::string(config.generatorSeqLen, chars[1]);
		}

		db.strings.push_back(seq1);
		db.seqs.push_back(db.strings.back());
		db.strings.push_back(seq2);
		db.seqs.push_back(db.strings.back());

		Comparison c{
			.originalComparisonIndex = i,
			.indexA = 2 * i,
			.sizeA = (int) db.seqs[2*i].size(),
			.indexB = 2 * i + 1,
			.sizeB = (int) db.seqs[2*i + 1].size(),
			.seeds = {{{seed_start, seed_start}, {seed_start, seed_start}}},
			.complexity = 0,
		};
		add_comparison<C>(db.cmps, c, swconfig.seedLength);
	}
	return std::move(db);
}

template SequenceDatabase<MultiComparison> loadSequences<MultiComparison>(GeneratorConfig& config, SWConfig& swconfig);
template SequenceDatabase<Comparison> loadSequences<Comparison>(GeneratorConfig& config, SWConfig& swconfig);
}