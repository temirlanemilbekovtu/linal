#include "linal/dim.hpp"

Dim::Dim(const size_t rows, const size_t cols)
    : rows(rows)
    , cols(cols) { }

Dim &Dim::operator=(const Dim &other) {
    if (this != &other) {
        rows = other.rows;
        cols = other.cols;
    }

    return *this;
}

size_t Dim::get_card() const {
    return rows * cols;
}
