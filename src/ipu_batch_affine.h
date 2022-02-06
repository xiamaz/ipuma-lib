#ifndef IPU_BATCH_AFFINE_HPP
#define IPU_BATCH_AFFINE_HPP

#include<vector>
#include "ipu_base.h"
#include "partition.h"

using namespace poplar;

namespace ipu {
namespace batchaffine {

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
  int getSequenceBufferSize8b() const;
  int getSequenceBufferSize32b() const;
  int getMetaBufferSize32b() const;
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
  static size_t getMetaOffset(const IPUAlgoConfig& config);

 public:
  IPUAlgoConfig algoconfig;

  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig);
  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id);

  std::string printTensors();

  static std::vector<std::tuple<int, int>> fillBuckets(const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B, int& err);
  std::vector<std::tuple<int, int>> fillBuckets(const std::vector<std::string>& A, const std::vector<std::string>& B, int& err);
  static void checkSequenceSizes(const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B);

  BlockAlignmentResults get_result();

  // Local Buffers
  void compare_local(const std::vector<std::string>& A, const std::vector<std::string>& B, bool errcheck = false);

  void refetch();

  // Remote bufffer
  void prepared_remote_compare(int32_t* inputs_begin,  int32_t* inputs_end, int32_t* results_begin, int32_t* results_end);
  static void prepare_remote(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B,  int32_t* inputs_begin,  int32_t* inputs_end, std::vector<int>& deviceMapping);
};
}  // namespace batchaffine
}  // namespace ipu
#endif  // IPU_BATCH_AFFINE_HPP