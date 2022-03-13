#ifndef IPU_CONFIG_H
#define IPU_CONFIG_H
#include "swatlib/swatlib.h"
#include "types.h"
#include "nlohmann/json.hpp"

namespace ipu {

using json = nlohmann::json;

struct SWConfig {
        int gapInit = 0;
        int gapExtend = -1;
        int matchValue = 1;
        int mismatchValue = -1;
        int ambiguityValue = -1;
        swatlib::Similarity similarity = swatlib::Similarity::nucleicAcid;
        swatlib::DataType datatype = swatlib::DataType::nucleicAcid;
};

struct IPUAlgoConfig {
  int tilesUsed = 1; // number of active vertices
  int maxAB = 300; // maximum length of a single comparison
  int maxBatches = 20; // maximum number of comparisons in a single batch
  int bufsize = 3000; // total size of buffer for A and B individually
  VertexType vtype = VertexType::cpp;
  Algorithm fillAlgo = Algorithm::fillFirst;
  bool forwardOnly = false; // do not calculate the start position of a match, this should approx 2x performance, as no reverse pass is needed
  bool useRemoteBuffer = false;
  int transmissionPrograms = 1;  // number of separate transmission programs, use only with remote!
  int ioTiles = 0;

  // this is maxbatches * num_vertices
  int getTotalNumberOfComparisons() const;

  // size of the sequence buffer with 32bits
  int getBufsize32b() const;
  int getTotalBufsize32b() const;
  int getMetaBufferSize32b() const;
  int getLenBufferSize32b() const;
  int getInputBufferSize32b() const;
};

void to_json(json& j, const SWConfig& c);

void from_json(const json& j, SWConfig& c);

void from_json(const json& j, Algorithm& a);

void to_json(json& j, const Algorithm& a);

void from_json(const json& j, VertexType& t);

void to_json(json& j, const VertexType& t);

void to_json(json& j, const IPUAlgoConfig& c);

void from_json(const json& j, IPUAlgoConfig& c);

}

#endif