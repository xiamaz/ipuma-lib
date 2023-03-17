#pragma once

#include <string_view>
#include "ipu_config.h"
#include "core/xdrop.hpp"


namespace cpu {

  std::tuple<std::vector<uint8_t>, std::vector<uint8_t>> createSeedPartsIPU(const std::string_view& sequence, int seedStart, int seedLength, swatlib::Encoding& enc) {
    auto seqRight = sequence.substr(seedStart + seedLength, sequence.size() - (seedStart + seedLength));
    auto seqLeft = sequence.substr(0, seedStart);
    std::string strRight(seqRight);
    std::string strLeft(seqLeft.rbegin(), seqLeft.rend());

    return {
      enc.encode(strLeft),
      enc.encode(strRight),
    };
  }

  class IpumaCpuAligner {
    const ipu::SWConfig& config;
  public:
    IpumaCpuAligner(const ipu::SWConfig& config) : config(config) {
    }

    int align(const std::string_view& seqH, const std::string_view& seqV, int32_t seedHStart, int32_t seedVStart, int32_t seedLen) {
      auto enc = swatlib::getEncoder(config.datatype);
      auto encH = enc.encode(std::string{seqH});
      auto encV = enc.encode(std::string{seqV});
      auto sim = swatlib::selectMatrix(config.similarity, config.matchValue, config.mismatchValue, config.ambiguityValue);
      auto score = ipumacore::xdrop::seed_extend_restricted_cpu<10, 1>(
        encH, seedHStart, 
        encV, seedVStart,
        seedLen,
        sim, (int) 20000 * 0.45);
      return score;
    }
  };
}