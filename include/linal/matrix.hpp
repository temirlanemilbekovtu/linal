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
    const int               _rows;
    const int               _cols;
    std::unique_ptr<T[]>    _data;

public:
    Matrix(const int rows, const int cols);
    Matrix(const Matrix &matrix);

    Matrix<T> operator+(const Matrix<T> &matrix);
    Matrix<T> operator*(const Matrix<T> &matrix);
    Matrix<T> operator*(const T &value);

    int     get_cardinality();
    int     get_rows();
    int     get_cols();

    T       get(const int row, const int col);
    void    set(const int row, const int col, const T &value);
};


#endif //LINAL_MATRIX_HPP
