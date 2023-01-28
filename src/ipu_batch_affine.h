#ifndef IPU_BATCH_AFFINE_HPP
#define IPU_BATCH_AFFINE_HPP

#define WORK_QUEUE_SIZE 100'000

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

using JobId = long;
struct BlockAlignmentResults {
  std::vector<int32_t> scores;
  std::vector<int32_t> a_range_result;
  std::vector<int32_t> b_range_result;
};

struct Batch {
  std::vector<int32_t> inputs;
  std::vector<int32_t> results;
  std::vector<int32_t> origin_comparison_index;
  swatlib::TickTock batchTick;
  BlockAlignmentResults get_result();
  static void transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end);
  static void transferResults(int32_t* results_begin, int32_t* results_end, int* mapping_begin, int* mapping_end, int32_t* scores_begin, int32_t* scores_end, int32_t* arange_begin, int32_t* arange_end, int32_t* brange_begin, int32_t* brange_end, int numComparisons);
};

struct Job {
  Job(): tick({}) {};
  void join();

  JobId id;
  msd::channel<int>* done_signal;
  uint64_t h2dCycles;
  uint64_t innerCycles;

  swatlib::TickTock jobTick;
  Batch* batch;
};



/**
 * Workflow:
 * 1. Langer sequenz wird präperiert -> variable sized vector ab `Batch`es als output. (Sorting passiert hier)
 * 2. Submitted -> Batch wird zum Job,.. Job der für 1:1 Batch zuständig ist
 * 3. Job beendet, batch enthält informationen.
 *    -> wir können im Batch UPC++ pounter magie machen.
 * 
 * 
 * Seqs = {"ctgaa", "aact", "ttgg", "aact" }
 * Comps = [
 *            (0, 1), // index = 0
 *            (3, 4)  // index = 1
 *          ]
 * 
 */

class SWAlgorithm : public IPUAlgorithm {
 private:
  msd::channel<Job*> work_queue;
  std::mutex tableMutex;
  std::map<JobId, Job*> resultTable;
  std::vector<std::thread> executor_procs;

 public:
  IPUAlgoConfig algoconfig;

  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig);
  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id, size_t ipuCount = 1, bool runExecutor = true);

  static void checkSequenceSizes(const IPUAlgoConfig& algoconfig, const std::vector<int>& SeqSizes);
  static void checkSequenceSizes(const IPUAlgoConfig& algoconfig, const RawSequences& Seqs);

  static void fillBuckets(Algorithm algo, partition::BucketMap& map, const RawSequences& A, const RawSequences& B, int offset = 0);
  static void fillMNBuckets(Algorithm algo, partition::BucketMap& map, const RawSequences& Seqs, const Comparisons& Cmps, int offset = 0);

  static void fill_input_buffer(const partition::BucketMap& map, const swatlib::DataType dtype, const IPUAlgoConfig& algoconfig, const RawSequences& Seqs, const Comparisons& Cmps, int32_t* inputs_begin, int32_t* inputs_end, int32_t* mapping);
  static void fill_input_buffer(const partition::BucketMap& map, const swatlib::DataType dtype, const IPUAlgoConfig& algoconfig, const RawSequences& A, const RawSequences& B, int32_t* inputs_begin, int32_t* inputs_end, int32_t* mapping);

  // Local Buffers
  // void compare_local(const std::vector<std::string>& A, const std::vector<std::string>& B, bool errcheck = false);
  // void compare_mn_local(const std::vector<std::string>& Seqs, const Comparisons& Cmps, bool errcheck = false);

  // Convenience (checkSequenceSizes, fillMNBuckets, fill_input_buffer)
  std::vector<Batch> create_batches(
    const SWConfig& swconfig, 
    const IPUAlgoConfig& algoconfig, 
    const std::vector<std::string>& Seqs,
    const Comparisons& Cmps
  );

  Job* async_submit(Batch* batch);
  void blocking_join(Job& job);

  void run_executor();
  ~SWAlgorithm();
};
}  // namespace batchaffine
}  // namespace ipu
#endif  // IPU_BATCH_AFFINE_HPP