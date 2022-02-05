// Based on: https://stackoverflow.com/a/31392034
#ifndef SWATLIB_TIMING_H
#define SWATLIB_TIMING_H

#include <cmath>
#include <string>
#include <chrono>
#include <algorithm>
#include <iostream>

namespace swatlib {

template<typename Duration = std::chrono::milliseconds,
         typename F,
         typename ... Args>
inline typename Duration::rep profile(F&& fun,  Args&&... args) {
    const auto beg = std::chrono::high_resolution_clock::now();
    std::forward<F>(fun)(std::forward<Args>(args)...);
    const auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<Duration>(end - beg).count();
}

/**
 * No affiliation with popular social media platforms.
 * 
 * Usage: First call tick and then call tock. Can be used multiple times simply by calling tick again.
 */
class TickTock {
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::microseconds accumulate_time = std::chrono::microseconds::zero();
    bool started = false;
    bool ended = false;
public:
    void tick() {
        start = std::chrono::high_resolution_clock::now();
        started = true;
        ended = false;
    }

    void tock() {
        if (!started)
            throw std::runtime_error("Not started.");
        end = std::chrono::high_resolution_clock::now();
        ended = true;

        accumulate_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    }

    int64_t accumulate_microseconds() {
        return accumulate_time.count();
    }

    template<typename Unit = std::chrono::milliseconds>
    typename Unit::rep duration() {
        if (!(started && ended))
            throw std::runtime_error("Not started and ended yet.");
        return std::chrono::duration_cast<Unit>(end - start).count();
    }
};

}

#endif // TIMING_H