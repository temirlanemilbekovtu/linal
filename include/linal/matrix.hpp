#pragma once

#ifndef LINAL_MATRIX_HPP
#define LINAL_MATRIX_HPP

#include <memory>
#include <stdexcept>
#include <iostream>

#include "export.hpp"
#include "my_type_traits.hpp"

template <typename T>
struct LINAL_API Matrix {
    static_assert(is_addable<T>::value, "Type must support addition");
    static_assert(is_multipliable<T>::value, "Type must support multiplication");

private:
    const int                   _rows;
    const int                   _cols;
    std::unique_ptr<T[]>        _data;

public:
    explicit Matrix(const int rows, const int cols) : _rows(rows), _cols(cols) {
        _data = std::make_unique<T[]>(rows * cols);
    }

    Matrix(const Matrix<T> &matrix) : Matrix (matrix._rows, matrix._cols) {
        int cardinality = get_cardinality();
        for(int i = 0; i < cardinality; ++i) {
            _data[i] = matrix._data[i];
        }
    }

    Matrix<T> operator+(const Matrix<T> &matrix) const {
        if (_rows != matrix._rows || _cols != matrix._cols) {
            throw std::invalid_argument("Matrix sizes don't match for the sum operation!");
        }

        Matrix<T> result(_rows, _cols);

        for(int i = 0; i < _rows; ++i) {
            for(int j = 0; j < _cols; ++j) {
                T entry = get_entry(i, j) + matrix.get_entry(i, j);
                result.set_entry(i, j, entry);
            }
        }

        return result;
    }

    Matrix<T> operator*(const Matrix<T> &matrix) const {
        if (_cols != matrix._rows) {
            throw std::invalid_argument("Matrix sizes don't match for the multiplication operation!");
        }

        return mul_classic(this, matrix);
    }

    Matrix<T> operator*(const T &value) const {
        Matrix<T> result(_rows, _cols);

        for(int i = 0; i < _rows; ++i) {
            for(int j = 0; j < _cols; ++j) {
                T entry = get_entry(i, j) * value;
                result.set_entry(i, j, entry);
            }
        }

        return result;
    }

    static Matrix<T> mul_classic(const Matrix<T> &a, const Matrix<T> &b) {
        if (a.get_cols() != b.get_rows()) {
            throw std::invalid_argument("Matrix sizes don't match for the multiplication operation!");
        }

        const int &m = a.get_rows();
        const int &n = a.get_cols();
        const int &k = b.get_cols();
        Matrix<T> result(m, k);

        for(int i = 0; i < k; ++i) {
            for(int j = 0; j < m; ++j) {
                T entry {};
                for(int l = 0; l < n; ++l) {
                    entry += a.get_entry(m, n) * b.get_entry(n, k);
                }
                result.set_entry(m, k, entry);
            }
        }

        return result;
    }

    static Matrix<T> mul_winograd(const Matrix<T> &a, const Matrix<T> &b) {
        if (a.get_cols() != b.get_rows()) {
            throw std::invalid_argument("Matrix sizes don't match for the multiplication operation!");
        }

        const int &m = a.get_rows();
        const int &n = a.get_cols();
        const int &k = b.get_cols();
        Matrix<T> result{m, k};

        int n_half = n / 2;
        int n_less = n - 1;
        T row_factor[m];
        T col_factor[k];

        for(int i = 0; i < m; ++i) {
            row_factor[i] = T();
            for(int j = 0; j < n_half; j += 2) {
                row_factor[i] += a.get_entry(i, j) * a.get_entry(i, j + 1);
            }
        }

        for(int i = 0; i < k; ++i) {
            col_factor[i] = T();
            for(int j = 0; j < n_half; j += 2) {
                col_factor[i] += b.get_entry(j, i) * b.get_entry(j + 1, i);
            }
        }

        for(int i = 0; i < m; ++i) {
            for(int j = 0; j < k; ++j) {
                T entry = -1 * (row_factor[i] + col_factor[j]);
                for(int t = 0; t < n_half; t += 2) {
                    entry += (a.get_entry(i, t) + b.get_entry(t + 1, j)) * (a.get_entry(i, t + 1) + b.get_entry(t, j));
                }
                result.set_entry(i, j, entry);
            }
        }

        if (n % 2 == 1) {
            for(int i = 0; i < m; ++i) {
                for(int j = 0; j < k; ++j) {
                    T entry = result.get_entry(i, j) + get_entry(i, n_less) * matrix.get_entry(n_less, j);
                    result.set_entry(i, j, entry);
                }
            }
        }

        return result;
    }

    int get_cardinality() const {
        return _rows * _cols;
    }

    int get_rows() const {
        return _rows;
    }

    int get_cols() const {
        return _cols;
    }

    T get_entry(const int row, const int col) const {
        int index = row * _cols + col;
        return _data[index];
    }

    void set_entry(const int row, const int col, const T &value) {
        int index = row * _cols + col;
        _data[index] = value;
    }
};

#endif //LINAL_MATRIX_HPP