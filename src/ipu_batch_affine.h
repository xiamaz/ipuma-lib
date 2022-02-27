#ifndef IPU_BATCH_AFFINE_HPP
#define IPU_BATCH_AFFINE_HPP

#include <vector>
#include "ipu_base.h"
#include "partition.h"
#include "ipu_config.h"


using namespace poplar;

namespace ipu {
namespace batchaffine {
struct BlockAlignmentResults {
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
};

using slotToken = int;
class SWAlgorithm : public IPUAlgorithm {
 private:
  // std::vector<int32_t> results;
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
  int thread_id;

  int slot_avail;
  int last_slot;
  std::vector<bool> slots;

  static size_t getSeqsOffset(const IPUAlgoConfig& config);
  static size_t getMetaOffset(const IPUAlgoConfig& config);


  void release_slot(slotToken i);

 public:
  IPUAlgoConfig algoconfig;
  bool use_remote_buffer;

  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig);
  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id, bool useRemoteBuffer = true, size_t slotCap = 1);

  std::string printTensors();

  bool slot_available();
  slotToken queue_slot();

  static void checkSequenceSizes(const IPUAlgoConfig& algoconfig, const std::vector<int>& SeqSizes);
  static void checkSequenceSizes(const IPUAlgoConfig& algoconfig, const RawSequences& Seqs);

  static void fillBuckets(Algorithm algo, partition::BucketMap& map, const RawSequences& A, const RawSequences& B, int offset = 0);
  static void fillMNBuckets(Algorithm algo, partition::BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int offset = 0);

  static void fill_input_buffer(const partition::BucketMap& map, const swatlib::DataType dtype, const IPUAlgoConfig& algoconfig, const RawSequences& Seqs, const Comparisons& Cmps, int32_t* inputs_begin, int32_t* inputs_end, int32_t* mapping);
  static void fill_input_buffer(const partition::BucketMap& map, const swatlib::DataType dtype, const IPUAlgoConfig& algoconfig, const RawSequences& A, const RawSequences& B, int32_t* inputs_begin, int32_t* inputs_end, int32_t* mapping);

  BlockAlignmentResults get_result();

  // Local Buffers
  void compare_local(const std::vector<std::string>& A, const std::vector<std::string>& B, bool errcheck = false);
  void compare_mn_local(const std::vector<std::string>& Seqs, const Comparisons& Cmps, bool errcheck = false);

  void refetch();

  // Remote bufffer
  void prepared_remote_compare(int32_t* inputs_begin,  int32_t* inputs_end, int32_t* results_begin, int32_t* results_end, slotToken slot_token = 0);

  void upload(int32_t* inputs_begin, int32_t* inputs_end, slotToken slot);
  slotToken upload(int32_t* inputs_begin, int32_t* inputs_end);
 
  static void prepare_remote(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B,  int32_t* inputs_begin, int32_t* inputs_end, int* deviceMapping);
  static void transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end);
  static void transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end, int numComparisons);
};
}  // namespace batchaffine
}  // namespace ipu
#endif  // IPU_BATCH_AFFINE_HPP