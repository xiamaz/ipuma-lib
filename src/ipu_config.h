#ifndef IPU_CONFIG_H
#define IPU_CONFIG_H
#include "swatlib/swatlib.h"
#include "types.h"
#include "nlohmann/json.hpp"

#ifndef KLIGN_IPU_MAX_SEQUENCE_LENGTH
#define KLIGN_IPU_MAX_SEQUENCE_LENGTH 10000
#endif

#ifndef KLIGN_IPU_COMPLEXITY_ALGO
#define KLIGN_IPU_COMPLEXITY_ALGO "xdrop"
#endif

#ifndef KLIGN_IPU_NUM_VERTICES
#define KLIGN_IPU_NUM_VERTICES 1472
#endif

#ifndef KLIGN_IPU_MAX_CMPS_PER_VERTEX
#define KLIGN_IPU_MAX_CMPS_PER_VERTEX 200
#endif

#ifndef KLIGN_IPU_VERTEX_BUFFER_SIZE
#define KLIGN_IPU_VERTEX_BUFFER_SIZE 120000
#endif

#ifndef KLIGN_IPU_VTYPE
#define KLIGN_IPU_VTYPE "multixdrop"
#endif

#ifndef KLIGN_IPU_PARTITION_ALGO
#define KLIGN_IPU_PARTITION_ALGO "greedy"
#endif

#ifndef ALN_GAP_OPENING_COST
#define ALN_GAP_OPENING_COST 1
#endif

#ifndef ALN_GAP_EXTENDING_COST
#define ALN_GAP_EXTENDING_COST 1
#endif

#ifndef ALN_MATCH_SCORE
#define ALN_MATCH_SCORE 1
#endif

#ifndef ALN_MISMATCH_COST
#define ALN_MISMATCH_COST 1
#endif

#ifndef ALN_AMBIGUITY_COST
#define ALN_AMBIGUITY_COST 1
#endif

#ifndef ALN_SIMILARITY
#define ALN_SIMILARITY "na"
#endif

#ifndef ALN_DATATYPE
#define ALN_DATATYPE "na"
#endif

namespace ipu {

using json = nlohmann::json;

struct SWConfig {
        int gapInit = -(ALN_GAP_OPENING_COST - ALN_GAP_EXTENDING_COST);
        int gapExtend = -ALN_GAP_EXTENDING_COST;
        int matchValue = ALN_MATCH_SCORE;
        int mismatchValue = -ALN_MISMATCH_COST;
        int ambiguityValue = -ALN_AMBIGUITY_COST;
        swatlib::Similarity similarity = swatlib::strToSimilarity(ALN_SIMILARITY);
        swatlib::DataType datatype = swatlib::strToDataType(ALN_DATATYPE);
};

struct IPUAlgoConfig {
  int numVertices = KLIGN_IPU_NUM_VERTICES; // number of active vertices
  int maxSequenceLength = KLIGN_IPU_MAX_SEQUENCE_LENGTH; // maximum length of a single comparison
  int maxComparisonsPerVertex = KLIGN_IPU_MAX_CMPS_PER_VERTEX; // maximum number of comparisons in a single batch
  int vertexBufferSize = KLIGN_IPU_VERTEX_BUFFER_SIZE; // total size of buffer for A and B individually
  VertexType vtype = strToVertexType(KLIGN_IPU_VTYPE);
  Algorithm fillAlgo = strToAlgorithm(KLIGN_IPU_PARTITION_ALGO);
  Complexity complexityAlgo = strToComplexity(KLIGN_IPU_COMPLEXITY_ALGO);
  bool partitioningSortComparisons = true;
  bool forwardOnly = false; // do not calculate the start position of a match, this should approx 2x performance, as no reverse pass is needed
  int ioTiles = 0;

  // Optional: XDrop
  int xDrop = 10;
  double bandPercentageXDrop = 0.5;
  int seedLength = -1;

  // this is maxbatches * num_vertices and is the maximum number of comparisons in a single batch
  int getTotalNumberOfComparisons() const;

  // size of the sequence buffer with 32bits
  int getBufsize32b() const;
  int getVertexBufsize32b() const;
  int getMetaBufferSize32b() const;
  int getLenBufferSize32b() const;
  int getInputBufferSize32b() const;

  int getMetaStructSize32b() const;

  size_t getOffsetInputSequence() const;
  size_t getOffsetMetadata() const;
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