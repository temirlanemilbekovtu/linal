#pragma once

#ifndef LINAL_MATRIX_HPP
#define LINAL_MATRIX_HPP

#include <memory>
#include "export.hpp"
#include "my_type_traits.hpp"

template <typename T>
struct LINAL_API Matrix {
    static_assert(is_addable<T>::value_type, "Type must support addition");
    static_assert(is_multipliable<T>::value_type, "Type must support multiplication");

private:
    const int                   _rows;
    const int                   _cols;
    std::unique_ptr<T[]>        _data;

public:
    explicit        Matrix(const int rows, const int cols);
                    Matrix(const Matrix<T> &matrix);

    Matrix<T>       operator+(const Matrix<T> &matrix) const;
    Matrix<T>       operator*(const Matrix<T> &matrix) const;
    Matrix<T>       operator*(const T &value) const;

    int             get_cardinality() const;
    int             get_rows() const;
    int             get_cols() const;

    T               get_entry(const int row, const int col) const;
    void            set_entry(const int row, const int col, const T &value);
};


#endif //LINAL_MATRIX_HPP
