#ifndef SLICE_H
#define SLICE_H

#include "pch.h"

class Slice {
public:
    size_t start, end, step;

    Slice(size_t start = 0, size_t end = 0, size_t step = 1)
        : start(start), end(end), step(step) {}

    // Check slicing range
    void check(size_t max_size) const {
        if (start >= max_size || end > max_size || step == 0)
            throw std::out_of_range("Invalid slice range.");
    }
};

#endif
