#ifndef IPU_MULTI_DRIVER_HPP
#define IPU_MULTI_DRIVER_HPP
#include "ipuswconfig.hpp"
#include "nlohmann/json.hpp"
#include <plog/Log.h>

using json = nlohmann::json;

class IPUMultiDriver {
        ipu::batchaffine::SWAlgorithm* driver;

        double clockFrequency;
public:
        IpuSwConfig config;
				std::atomic<int> totalCmpsProcessed = 0;
				std::atomic<int> totalBatchesProcessed = 0;

        IPUMultiDriver(IpuSwConfig config) : config(config) {
                PLOGD << "Initialize drivers.";
                driver = new ipu::batchaffine::SWAlgorithm(config.swconfig, config.ipuconfig, 0, 0, config.numDevices);
                clockFrequency = driver->getTileClockFrequency();
        }

        ~IPUMultiDriver() {
                delete driver;
        }

        void fill_input_buffer(const ipu::RawSequences& a, const ipu::RawSequences& b, std::vector<int32_t>& input_buffer, std::vector<int32_t>& mapping) {
                ipu::batchaffine::SWAlgorithm::prepare_remote(config.swconfig, config.ipuconfig, a, b, &*input_buffer.begin(), &*input_buffer.end(), mapping.data());
        }

        void submitWait(std::vector<int32_t>& inputBuffer, std::vector<int32_t>& resultBuffer) {
                auto* job = driver->async_submit_prepared_remote_compare(&*inputBuffer.begin(), &*inputBuffer.end(), &*resultBuffer.begin(), &*resultBuffer.end());
                job->join();

        PLOGD << "Total engine run time (in s): " << static_cast<double>(job->tick.duration<std::chrono::milliseconds>()) / 1000.0;
#ifdef IPUMA_DEBUG
                // auto [lot, sid] = unpack_slot(slot_token);
                auto cyclesOuter = job->h2dCycles + job->innerCycles;
                auto cyclesInner = job->innerCycles;
								auto timeJob = static_cast<double>(job->tick.accumulate_microseconds()) / 1e6;
								auto timeBatch = static_cast<double>(job->sb.runTick.accumulate_microseconds()) / 1e6;
                auto timeOuter = static_cast<double>(cyclesOuter) / clockFrequency;
                auto timeInner = static_cast<double>(cyclesInner) / clockFrequency;

                int32_t *meta_input = job->sb.inputs_begin + config.ipuconfig.getTotalBufsize32b();

                // GCUPS computation
                uint64_t cellCount = 0;
                uint64_t dataCount = 0;
                for (size_t i = 0; i < config.ipuconfig.getTotalNumberOfComparisons(); i++)
                {
                        auto a_len = meta_input[4 * i];
                        auto b_len = meta_input[4 * i + 2];
                        cellCount += a_len * b_len;
                        dataCount += a_len + b_len;
                }

                double GCUPSOuter = static_cast<double>(cellCount) / timeOuter / 1e9;
                double GCUPSInner = static_cast<double>(cellCount) / timeInner / 1e9;
                double GCUPSJob = static_cast<double>(cellCount) / timeJob / 1e9;
                double GCUPSBatch = static_cast<double>(cellCount) / timeBatch / 1e9;

                double totalTransferSize = config.ipuconfig.getInputBufferSize32b() * 4;

                auto transferTime = timeOuter - timeInner;
                auto transferInfoRatio = static_cast<double>(dataCount) / totalTransferSize * 100;
                auto transferBandwidth = totalTransferSize / transferTime / 1e6;
                auto transferBandwidthPerVertex = transferBandwidth / config.ipuconfig.tilesUsed;

								json log = {
									{"tag", "batch_perf"},
									{"cycle_inner", cyclesInner},
									{"cycle_outer", cyclesOuter},
									{"time_inner", timeInner},
									{"time_outer", timeOuter},
									{"time_job", timeJob},
									{"cell_count", cellCount},
									{"gcups_inner", GCUPSInner},
									{"gcups_outer", GCUPSOuter},
									{"gcups_job", GCUPSJob},
									{"gcups_batch", GCUPSBatch},
									{"transfer_size_total", totalTransferSize},
									{"time_h2d", transferTime},
									{"data_count", dataCount},
									{"transfer_ratio", transferInfoRatio},
									{"transfer_bandwidth", transferBandwidth},
									{"transfer_bandwidth_per_vertex", transferBandwidthPerVertex},
								};
								PLOGD << IPU_JSON_LOG_TAG << log.dump();
#endif
                delete job;
        }

        void fillResults(std::vector<int32_t>& resultBuffer, std::vector<int>& mappingBuffer, ipu::batchaffine::BlockAlignmentResults& results, const int numCmps) {
                results.scores.resize(numCmps);
                results.a_range_result.resize(numCmps);
                results.b_range_result.resize(numCmps);
    ipu::batchaffine::SWAlgorithm::transferResults(
        &*resultBuffer.begin(), &*resultBuffer.end(),
        &*mappingBuffer.begin(), &*mappingBuffer.end(),
        &*(results.scores.begin()), &*(results.scores.begin() + numCmps),
        &*(results.a_range_result.begin()), &*(results.a_range_result.begin() + numCmps),
        &*(results.b_range_result.begin()), &*(results.b_range_result.begin() + numCmps), numCmps);
        }
};

#endif