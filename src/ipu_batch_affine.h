#ifndef IPU_BATCH_AFFINE_HPP
#define IPU_BATCH_AFFINE_HPP

#include <vector>
#include <thread>
#include <map>


#include "ipu_base.h"
#include "ipu_config.h"
#include "partition.h"
#include "msd/channel.hpp"

using namespace poplar;

namespace ipu {
namespace batchaffine {

using slotToken = long;
struct BlockAlignmentResults {
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
};

struct SubmittedBatch {
  int32_t* inputs_begin;
  int32_t* inputs_end;
  int32_t* results_begin;
  int32_t* results_end;
  slotToken slot;
  msd::channel<int>* signal_done;
  uint64_t* cyclesH2D = nullptr;
  uint64_t* cyclesInner = nullptr;
  swatlib::TickTock runTick;
};

struct Job {
  Job(): tick({}) {};
  void join();
  msd::channel<int>* done_signal;
  uint64_t h2dCycles;
  uint64_t innerCycles;
  swatlib::TickTock tick;
  SubmittedBatch sb;
};

using BatchChannel = msd::channel<SubmittedBatch*>;

class SWAlgorithm : public IPUAlgorithm {
 private:
  msd::channel<SubmittedBatch*> work_queue;
  std::mutex tableMutex;
  std::map<slotToken, SubmittedBatch*> resultTable;
  // std::vector<int32_t> results;
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
  std::vector<std::thread> executor_procs;

  int slot_size;
  std::vector<int> slot_avail;
  int last_slot;
  std::vector<std::vector<bool>> slots;

  static size_t getSeqsOffset(const IPUAlgoConfig& config);
  static size_t getMetaOffset(const IPUAlgoConfig& config);

  void release_slot(slotToken i);

 public:
  IPUAlgoConfig algoconfig;

  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig);
  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id, size_t slotCap = 1, size_t ipuCount = 1);

  std::string printTensors();

  static int calculate_slot_region_index(const IPUAlgoConfig& algoconfig, int max_buffer_size);
  int calculate_slot_region_index(int max_buffer_size);
  std::tuple<int, slotToken> unpack_slot(slotToken);
  bool slot_available(int max_buffer_size);
  slotToken queue_slot(int max_buffer_size);

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
  void prepared_remote_compare(int32_t* inputs_begin, int32_t* inputs_end, int32_t* results_begin, int32_t* results_end, slotToken slot_token = 0);

  Job* async_submit_prepared_remote_compare(int32_t* inputs_begin, int32_t* inputs_end, int32_t* results_begin, int32_t* results_end);
  void blocking_join_prepared_remote_compare(Job& job);

  void upload(int32_t* inputs_begin, int32_t* inputs_end, slotToken slot);
  // slotToken upload(int32_t* inputs_begin, int32_t* inputs_end);

  void run_executor();
  static int prepare_remote(const SWConfig& swconfig, const IPUAlgoConfig& algoconfig, const std::vector<std::string>& A, const std::vector<std::string>& B, int32_t* inputs_begin, int32_t* inputs_end, int* deviceMapping);
  static void transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end);
  static void transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end, int numComparisons);
  ~SWAlgorithm();
};
}  // namespace batchaffine
}  // namespace ipu
#endif  // IPU_BATCH_AFFINE_HPP