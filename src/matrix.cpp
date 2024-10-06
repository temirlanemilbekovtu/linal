#include <stdexcept>

#include "linal/matrix.hpp"

template<typename T>
Matrix<T>::Matrix(const int rows, const int cols)
: _rows(rows), _cols(cols) {
    _data = std::make_unique<T[]>(rows * cols);
}

template<typename T>
Matrix<T>::Matrix(const Matrix &matrix)
: _rows(matrix._rows), _cols(matrix._cols) {
    int cardinality = get_cardinality();
    for(int i = 0; i < cardinality; ++i) {
        _data[i] = matrix._data[i];
    }
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &matrix) {
    if (_rows != matrix._rows || _cols != matrix._cols) {
        throw std::invalid_argument("Matrix sizes don't match for the sum operation!");
    }

    Matrix<T> result(_rows, _cols);

    for(int i = 0; i < _rows; ++i) {
        for(int j = 0; j < _cols; ++j) {
            T entry = get(i, j) + matrix.get(i, j);
            result.set(i, j, entry);
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &matrix) {
    if (_cols != matrix._rows) {
        throw std::invalid_argument("Matrix sizes don't match for the multiplication operation!");
    }

    const int &m = _rows;
    const int &n = _cols;
    const int &k = matrix._cols;
    Matrix<T> result(m, k);

    auto compute_factor = [](Matrix<T> &matrix, T *factor[], int base, int n) {
        for(int i = 0; i < base; ++i) {
            factor[i] = 0;
            for(int j = 0; j < n; j += 2) {
                 factor += matrix.get(i, j) * matrix.get(i, j + 1);
            }
        }
    };

    T row_factor[m];
    T col_factor[k];

    compute_factor(this, row_factor, m, n);
    compute_factor(matrix, col_factor, k, n);

    for(int i = 0; i < m; ++i) {
        for(int j = 0; j < k; ++j) {
            T entry = (row_factor[i] * col_factor[j]) * -1;
            for(int t = 0; t < n; t += 2) {
                entry += (get(i, t) + matrix.get(t + 1, j)) * (get(i, t + 1) + matrix.get(t, j));
            }
            result.set(i, j, entry);
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T &value) {
    Matrix<T> result(_rows, _cols);

    for(int i = 0; i < _rows; ++i) {
        for(int j = 0; j < _cols; ++j) {
            T entry = get(i, j) * value;
            result.set(i, j, entry);
        }
    }

    return result;
}


template<typename T>
int Matrix<T>::get_cardinality() {
    return _rows * _cols;
}

template<typename T>
int Matrix<T>::get_rows() {
    return _rows;
}

template<typename T>
int Matrix<T>::get_cols() {
    return _cols;
}

template<typename T>
T Matrix<T>::get(const int row, const int col) {
    int index = row * _cols + col;
    return _data[index];
}

template<typename T>
void Matrix<T>::set(const int row, const int col, const T &value) {
    int index = row * _cols + col;
    _data[index] = value;
}

