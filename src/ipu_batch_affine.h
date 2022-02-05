#ifndef IPU_BATCH_AFFINE_HPP
#define IPU_BATCH_AFFINE_HPP

#include "ipu_base.h"
#include<vector>

using namespace poplar;

namespace ipu {
namespace batchaffine {

namespace partition {
  enum class Algorithm {fillFirst, roundRobin, greedy};

  struct BucketData {
    int count;
    int lenA;
    int lenB;
    int weight;
  };

  int fillFirst(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity);
  int roundRobin(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity);
  int greedy(std::vector<std::tuple<int, int>>& mapping, const std::vector<std::string>& A, const std::vector<std::string>& B, int bucketCount, int bucketCapacity, int bucketCountCapacity);
}

const std::string STREAM_A = "a-write";
const std::string STREAM_A_LEN = "a-len-write";
const std::string STREAM_B = "b-write";
const std::string STREAM_B_LEN = "b-len-write";
const std::string STREAM_SCORES = "scores-read";
const std::string STREAM_MISMATCHES = "mismatches-read";
const std::string STREAM_A_RANGE = "a-range-read";
const std::string STREAM_B_RANGE = "b-range-read";
const std::string STREAM_CONCAT_ALL = "concat-read-all";
const std::string HOST_STREAM_CONCAT = "host-stream-concat";

const std::string IPU_AFFINE_CPP = "SWAffine";
const std::string IPU_AFFINE_ASM = "SWAffineAsm";

enum class VertexType { cpp, assembly, multi, multiasm, stripedasm, multistriped, multistripedasm };

static const std::string typeString[] = {"SWAffine", "SWAffineAsm", "MultiSWAffine", "MultiSWAffineAsm", "StripedSWAffineAsm", "MultiSWStriped", "MultiSWStripedAsm"};
std::string vertexTypeToString(VertexType v);

struct IPUAlgoConfig {
  int tilesUsed = 1; // number of active vertices
  int maxAB = 300; // maximum length of a single comparison
  int maxBatches = 20; // maximum number of comparisons in a single batch
  int bufsize = 3000; // total size of buffer for A and B individually
  VertexType vtype = VertexType::cpp;
  partition::Algorithm fillAlgo = partition::Algorithm::fillFirst;

  /**
   * @brief This calculates the total number of comparisons that can be computed in a single engine run on the IPU.
   * 
   * @return int 
   */
  int getTotalNumberOfComparisons() const;

  /**
   * @brief This calculated the required buffer size for input strings across all vertices.
   * 
   * @return int 
   */
  int getTotalBufferSize8b() const;
  int getTotalBufferSize32b() const;
  int getLenBufferSize32b() const;
  int getInputBufferSize32b() const;
};

struct BlockAlignmentResults {
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
};

class SWAlgorithm : public IPUAlgorithm {
 private:
  // std::vector<int32_t> results;
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
  int thread_id;

  static size_t getAOffset(const IPUAlgoConfig& config);
  static size_t getBOffset(const IPUAlgoConfig& config);
  static size_t getAlenOffset(const IPUAlgoConfig& config);
  static size_t getBlenOffset(const IPUAlgoConfig& config);

 public:
  IPUAlgoConfig algoconfig;

  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig);
  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id);

  std::string printTensors();

  static std::vector<std::tuple<int, int>> fillBuckets(IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B, int& err);
  std::vector<std::tuple<int, int>> fillBuckets(const std::vector<std::string>& A, const std::vector<std::string>& B, int& err);
  static void checkSequenceSizes(IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B);

  BlockAlignmentResults get_result();

  // Local Buffers
  void compare_local(const std::vector<std::string>& A, const std::vector<std::string>& B, bool errcheck = true);

  void refetch();

  // Remote bufffer
  void prepared_remote_compare(int32_t* inputs_begin,  int32_t* inputs_end, int32_t* results_begin, int32_t* results_end);
  static void prepare_remote(IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B,  int32_t* inputs_begin,  int32_t* inputs_end, std::vector<int>& deviceMapping);
};
}  // namespace batchaffine
}  // namespace ipu
#endif  // IPU_BATCH_AFFINE_HPP