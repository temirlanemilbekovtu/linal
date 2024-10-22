#pragma once

#ifndef LINAL_DIM_HPP
#define LINAL_DIM_HPP

#include <cstddef>

#include "export.hpp"

struct LINAL_API Dim {
    size_t rows;
    size_t cols;

    explicit Dim(const size_t rows, const size_t cols);

    Dim &operator=(const Dim &other);
    bool operator==(const Dim &other) const;
    bool operator!=(const Dim &other) const;

    size_t get_card() const;
};

#endif //LINAL_DIM_HPP
