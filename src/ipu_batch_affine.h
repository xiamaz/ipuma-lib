#ifndef IPU_BATCH_AFFINE_HPP
#define IPU_BATCH_AFFINE_HPP

#define WORK_QUEUE_SIZE 100'000

#include <vector>
#include <thread>
#include <map>


#include "ipu_base.h"
#include "ipu_config.h"
#include "batching.h"
#include "partition.h"
#include "msd/channel.hpp"

using namespace poplar;

namespace ipu {
namespace batchaffine {

using JobId = long;

struct Job {
  Job(): tick({}) {};
  void join();

  JobId id;
  msd::channel<int>* done_signal;
  uint64_t h2dCycles;  // cycles necessary for transfer host to device
  uint64_t innerCycles; // cycles necessary for graph computation

  swatlib::TickTock tick;
  Batch* batch;
};
using JobMap = std::map<JobId, Job*>;

void computeJobMetrics(const Job&, double, size_t, size_t);

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
 protected:
  IPUAlgoConfig algoconfig;

  static void checkSequenceSizes(const IPUAlgoConfig& algoconfig, const RawSequences& Seqs);
 public:

  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig);
  SWAlgorithm(SWConfig config, IPUAlgoConfig algoconfig, int thread_id, size_t ipuCount = 1, bool runExecutor = true);

  // Convenience (checkSequenceSizes, fillMNBuckets, fill_input_buffer)
  std::vector<Batch> create_batches(
    const RawSequences& Seqs,
    Comparisons& Cmps
  );

  Job* async_submit(Batch* batch);
  void blocking_join(Job& job);

  void run_executor();
  ~SWAlgorithm();
};
}  // namespace batchaffine
}  // namespace ipu
#endif  // IPU_BATCH_AFFINE_HPP