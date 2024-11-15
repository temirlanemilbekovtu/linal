#pragma once

#ifndef LINAL_DIM_HPP
#define LINAL_DIM_HPP

#include <cstddef>

#include "export.hpp"

struct LINAL_API Dim {
private:
    size_t rows;
    size_t cols;
    size_t card;

    void update_card();

public:
    Dim();
    explicit Dim(const size_t rows, const size_t cols);

    Dim &operator=(const Dim &other);
    bool operator==(const Dim &other) const;
    bool operator!=(const Dim &other) const;

    size_t get_rows() const;
    void   set_rows(const size_t value);
    size_t get_cols() const;
    void   set_cols(const size_t value);
    size_t get_card() const;
};

#endif //LINAL_DIM_HPP
