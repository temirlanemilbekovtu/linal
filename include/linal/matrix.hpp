#pragma once

#ifndef LINAL_MATRIX_HPP
#define LINAL_MATRIX_HPP

#include <memory>
#include <stdexcept>

#include "export.hpp"
#include "my_type_traits.hpp"
#include "utils.hpp"
#include "point.hpp"
#include "dim.hpp"

template <typename T>
struct LINAL_API Matrix {
public:
    static_assert(is_addable<T>::value, "Type must support addition!");
    static_assert(is_deductible<T>::value, "Type must support deduction!");
    static_assert(is_multipliable<T>::value, "Type must support multiplication!");
    static_assert(std::is_default_constructible<T>::value, "Type must define default constructor!");

    enum MulType {
        AUTO,
        CLASSIC,
        WINOGRAD,
        STRASSEN,
    };

    enum Quadrant {
        TOP_LEFT,
        TOP_RIGHT,
        BOTTOM_LEFT,
        BOTTOM_RIGHT,
    };

private:
    const Dim _size;
    std::unique_ptr<T[]> _data;

    static Matrix<T> mul_classic(const Matrix<T> &lhs, const Matrix<T> &rhs);
    static Matrix<T> mul_winograd(const Matrix<T> &lhs, const Matrix<T> &rhs);
    static Matrix<T> mul_strassen(const Matrix<T> &lhs, const Matrix<T> &rhs);

public:

    explicit Matrix(const Dim size);
    explicit Matrix(const Dim size, const T &default_value);
    Matrix(const Matrix<T> &other);

    static Matrix<T> sum(const Matrix<T> &lhs, const Matrix<T> &rhs);
    static Matrix<T> ddt(const Matrix<T> &lhs, const Matrix<T> &rhs);
    static Matrix<T> mul(const Matrix<T> &lhs, const Matrix<T> &rhs, MulType type);
    static Matrix<T> mul(const Matrix<T> &obj, const T &value);

    static void check_same_size(const Matrix<T> &a, const Matrix<T> &b);
    static void check_multiplicity(const Matrix<T> &lhs, const Matrix<T> &rhs);

    void add(const Matrix<T> &rhs);
    void ddt(const Matrix<T> &rhs);
    void mul(const T &value);

    Matrix<T> &operator=(const Matrix<T> &source);
    Matrix<T> operator+(const Matrix<T> &other) const;
    Matrix<T> operator-(const Matrix<T> &other) const;
    Matrix<T> operator*(const Matrix<T> &other) const;
    Matrix<T> operator*(const T &value) const;
    T* operator[](const size_t index);
    const T* operator[](const size_t index) const;

    Matrix<T> &get_square_power_of_two() const;

    bool is_square() const;
    bool is_square_power_of_two() const;
    Dim get_size() const;

    T &at(const Point index);
    const T &at(const Point index) const;
};

#pragma region constructors

template<typename T>
Matrix<T>::Matrix(const Dim size)
: _size(size)
, _data(std::make_unique<T[]>(size.get_card())) { }

template<typename T>
Matrix<T>::Matrix(const Dim size, const T &default_value) : Matrix(size) {
    std::fill(_data, _data + size.get_card(), default_value);
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T> &other) : Matrix(other._size) {
    int card = _size.get_card();
    for(int i = 0; i < card; ++i) {
        this[i] = other[i];
    }
}

#pragma endregion

#pragma region operators

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &source) {
    if (this != &source) {
        _size = source.get_size();
        size_t card = _size.get_card();
        _data = std::make_unique<T[]>(card);
        for(int i = 0; i < card; ++i) {
            this[i] = source[i];
        }
    }

    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &other) const {
    return sum(*this, other);
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &other) const {
    return ddt(*this, other);
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &other) const {
    // multiplication should be here
    return Matrix{other.get_size()};
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T &value) const {
    return mul(*this, value);
}

#pragma endregion

#pragma region arithmetic_static

template<typename T>
Matrix<T> Matrix<T>::sum(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    check_same_size(lhs, rhs);

    Matrix<T> result(lhs.get_size());

    size_t card = lhs.get_size().get_card();
    for(int i = 0; i < card; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::ddt(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    check_same_size(lhs, rhs);

    Matrix<T> result(lhs.get_size());

    size_t card = lhs.get_size().get_card();
    for(int i = 0; i < card; ++i) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::mul(const Matrix<T> &lhs, const Matrix<T> &rhs, MulType type) {
    switch (type) {
        case MulType::AUTO:
            break;
        case MulType::CLASSIC:
            return mul_classic(lhs, rhs);
        case MulType::WINOGRAD:
            return mul_classic(lhs, rhs);
        case MulType::STRASSEN:
            return mul_classic(lhs, rhs);
    }

    return mul_classic(lhs, rhs);
}

template<typename T>
Matrix<T> Matrix<T>::mul(const Matrix<T> &obj, const T &value) {
    Dim size = obj.get_size();
    Matrix<T> result{size};
    for(int i = 0; i < size.rows; ++i) {
        for(int j = 0; j < size.cols; ++j) {
            result[i][j] = obj[i][j] * value;
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::mul_classic(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    check_multiplicity(lhs, rhs);

    const int &m = lhs.get_size().rows;
    const int &n = lhs.get_size().cols;
    const int &k = rhs.get_size().cols;
    Matrix<T> result{Dim(m, k)};

    for(int i = 0; i < m; ++i) {
        for(int j = 0; j < k; ++j) {
            for(int l = 0; l < n; ++l) {
                result[i][j] += lhs[i][l] * rhs[l][j];
            }
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::mul_winograd(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    check_multiplicity(lhs, rhs);

    const size_t &m = lhs.get_size().rows;
    const size_t &n = lhs.get_size().cols;
    const size_t &k = rhs.get_size().cols;
    Matrix<T> result{Dim(m, k)};

    size_t n_half = n / 2;
    size_t n_less = n - 1;
    T row_factor[m] = {};
    T col_factor[k] = {};

    for(int i = 0; i < n_half; ++i) {
        int i_plus = i + 1;
        for(int j = 0; j < m; ++j) {
            row_factor[j] += lhs[j][i] * lhs[j][i_plus];
        }
        for(int j = 0; j < k; ++j) {
            col_factor[j] += rhs[i][j] * rhs[i_plus][j];
        }
    }

    for(int i = 0; i < m; ++i) {
        for(int j = 0; j < k; ++j) {
            for(int t = 0; t < n_half; t += 2) {
                result[i][j] += (lhs[i][t] + rhs[t + 1][j]) * (lhs[i][t + 1] + rhs[t][j]);
            }
            result[i][j] -= row_factor[i] + col_factor[j];
        }
    }

    if (n % 2 == 1) {
        for(int i = 0; i < m; ++i) {
            for(int j = 0; j < k; ++j) {
                result[i][j] += lhs[i][n_less] * rhs[n_less][j];
            }
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::mul_strassen(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    // to be implemented
}

#pragma endregion

#pragma region other_static

template<typename T>
void Matrix<T>::check_same_size(const Matrix<T> &a, const Matrix<T> &b) {
    if (a.get_size() != a.get_size()) {
        throw std::invalid_argument("Matrices must have the same dimensions!");
    }
}

template<typename T>
void Matrix<T>::check_multiplicity(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    if (lhs.get_size().cols != rhs.get_size().rows) {
        throw std::invalid_argument("Dimensions don't match for multiplication!");
    }
}

#pragma endregion

#pragma region arithmetic

template<typename T>
void Matrix<T>::add(const Matrix<T> &rhs) {
    check_same_size(*this, rhs);

    size_t card = _size.get_card();
    for(int i = 0; i < card; ++i) {
        this[i] += rhs[i];
    }
}

template<typename T>
void Matrix<T>::ddt(const Matrix<T> &rhs) {
    check_same_size(*this, rhs);

    size_t card = _size.get_card();
    for(int i = 0; i < card; ++i) {
        this[i] -= rhs[i];
    }
}

template<typename T>
void Matrix<T>::mul(const T &value) {
    size_t card = _size.get_card();
    for(int i = 0; i < card; ++i) {
       this[i] *= value;
    }
}

#pragma endregion

#pragma region other

template<typename T>
bool Matrix<T>::is_square() const {
    return _size.rows == _size.cols;
}

template<typename T>
bool Matrix<T>::is_square_power_of_two() const {
    return _size.rows == _size.cols && Utils::is_power_of_two(_size.rows);
}

template<typename T>
Dim Matrix<T>::get_size() const {
    return _size;
}

template<typename T>
T &Matrix<T>::at(const Point index) {
    return &_data[index.row * _size.cols + index.col];
}

template<typename T>
const T &Matrix<T>::at(const Point index) const {
    return &_data[index.row * _size.cols + index.col];
}

#pragma endregion

#endif //LINAL_MATRIX_HPP