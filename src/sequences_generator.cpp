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

SequenceGenerator::SequenceGenerator(GeneratorConfig config, SWConfig swconfig) : config{config} {
  std::mt19937 gen(config.seed);

	int seed_start = config.sequenceLength / 2 - swconfig.seedLength / 2;
	if (seed_start + swconfig.seedLength >= config.sequenceLength) {
		throw std::runtime_error("Seed does not fit into sequence length");
	}
	for (int i = 0; i < config.sequenceCount; ++i) {
		auto seq1 = generateRandomSequence(config.sequenceLength, swconfig.datatype, gen);
		auto seq2 = mutateSequence(seq1, swconfig.datatype, gen, config.similarity, seed_start, swconfig.seedLength);

		sequences.push_back(seq1);
		seqs.push_back(sequences.back());
		sequences.push_back(seq2);
		seqs.push_back(sequences.back());

		cmps.push_back({
			.originalComparisonIndex = i,
			.indexA = 2 * i,
			.sizeA = (int) seqs[2*i].size(),
			.indexB = 2 * i + 1,
			.sizeB = (int) seqs[2*i + 1].size(),
			.seeds = {{{seed_start, swconfig.seedLength}, {seed_start, swconfig.seedLength}}},
			.complexity = 0,
		});
	}
}

std::tuple<RawSequences, ipu::Comparisons> SequenceGenerator::get() {
  return {std::move(seqs), std::move(cmps)};
}

void to_json(json& j, const GeneratorConfig& c) {
	j = json {
		{"generatorLength", c.sequenceLength},
		{"generatorCount", c.sequenceCount},
		{"generatorSimilarity", c.similarity},
		{"generatorSeed", c.seed},
	};
}
void from_json(const json& j, GeneratorConfig& c) {
	j.at("generatorLength").get_to(c.sequenceLength);
	j.at("generatorCount").get_to(c.sequenceCount);
	j.at("generatorSimilarity").get_to(c.similarity);
	j.at("generatorSeed").get_to(c.seed);
}
}