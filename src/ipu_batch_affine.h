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
  int getBufsize32b() const;
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

struct Comparison {
  int cmpIndex;
  int indexA;
  int indexB;
};

struct BucketMapping {
  std::vector<Comparison> comparisons;
};

using slotToken = int;
class SWAlgorithm : public IPUAlgorithm {
 private:
  // std::vector<int32_t> results;
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
  int thread_id;
  bool use_remote_buffer;

  int buf_size;
  int buf_cap;
  int last_slot;
  std::vector<bool> slots;

  static size_t getSeqsOffset(const IPUAlgoConfig& config);
  static size_t getMetaOffset(const IPUAlgoConfig& config);

  slotToken queue_slot() {
    assert(buf_has_capacity);
    int s = -1;
    for (size_t i = 0; i < slots.size(); i++) {
       if (!slots[i]) {
        s = i;
        break;
      }
    }
    assert(s != -1);
    slots[s] = true;
    buf_cap--;
    last_slot = s;
    return s;
  }

  void release_slot(slotToken i) {
    slots[i] = false;
    buf_cap++;
  }

 public:
  IPUAlgoConfig algoconfig;

  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig);
  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id, bool useRemoteBuffer = true, size_t bufSize = 1);

  std::string printTensors();

  bool buf_has_capacity() {
    return buf_cap > 0;
  }


  static std::vector<std::tuple<int, int>> fillBuckets(const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B, int& err);
  static std::vector<BucketMapping> fillMNBuckets(const IPUAlgoConfig& algoconfig, const std::vector<std::string>& Seqs, const std::vector<int>& comparisons);
  static void checkSequenceSizes(const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B);

  BlockAlignmentResults get_result();

  // Local Buffers
  void compare_local(const std::vector<std::string>& A, const std::vector<std::string>& B, bool errcheck = false);
  void compare_mn_local(const std::vector<std::string>& Seqs, const std::vector<int>& comparisons, bool errcheck = false);

  void refetch();

  // Remote bufffer
  void prepared_remote_compare(int32_t* inputs_begin,  int32_t* inputs_end, int32_t* results_begin, int32_t* results_end, slotToken slot_token = 0);
  slotToken upload(int32_t* inputs_begin, int32_t* inputs_end);
 
  static void prepare_remote(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B,  int32_t* inputs_begin, int32_t* inputs_end, int* deviceMapping);
  static std::vector<int> fill_input_buffer(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& Seqs, const std::vector<BucketMapping>& comparisonMapping, int numComparisons, int32_t* inputs_begin, int32_t* inputs_end);
  static void transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end);
};
}  // namespace batchaffine
}  // namespace ipu
#endif  // IPU_BATCH_AFFINE_HPP