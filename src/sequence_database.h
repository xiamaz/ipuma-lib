#pragma once
#include "types.h"

namespace ipu {

class SequenceDatabase {
public:
	virtual std::tuple<ipu::RawSequences, ipu::Comparisons> get() = 0;
};

}