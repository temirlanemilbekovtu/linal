#include "linal/dim.hpp"

Dim::Dim() : Dim(0, 0) { }

Dim::Dim(const size_t rows, const size_t cols)
    : rows(rows)
    , cols(cols)
    , card(rows * cols) { }

Dim &Dim::operator=(const Dim &other) {
    if (this != &other) {
        rows = other.rows;
        cols = other.cols;
    }

    return *this;
}

bool Dim::operator==(const Dim &other) const {
    return rows == other.rows && cols == other.cols;
}

bool Dim::operator!=(const Dim &other) const {
    return rows != other.rows || cols != other.cols;
}

void Dim::update_card() {
    card = rows * cols;
}

size_t Dim::get_rows() const {
    return rows;
}

void Dim::set_rows(const size_t value) {
    rows = value;
    update_card();
}

size_t Dim::get_cols() const {
    return cols;
}

void Dim::set_cols(const size_t value) {
    cols = value;
    update_card();
}

size_t Dim::get_card() const {
    return card;
}
