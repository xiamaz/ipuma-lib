#pragma once
#include "types.h"

namespace ipu {

template<typename C>
struct SequenceDatabase {
	std::vector<std::string> strings;
	ipu::RawSequences seqs;
	std::vector<C> cmps;

	std::tuple<ipu::RawSequences, std::vector<C>> get() {
  	return {std::move(seqs), std::move(cmps)};
	};
};

}