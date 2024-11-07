#pragma once

#ifndef LINAL_MATRIX_HPP
#define LINAL_MATRIX_HPP

#include <memory>
#include <array>
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
    static_assert(std::is_copy_assignable<T>::value, "Type must support copy assignment!");
    static_assert(std::is_default_constructible<T>::value, "Type must define default constructor!");

    enum MulType {
        AUTO,
        CLASSIC,
        WINOGRAD,
        STRASSEN,
    };

private:
    Dim _size;
    Point _offset;
    size_t _offset_plain;
    size_t _total_cols;
    std::shared_ptr<T[]> _data;

    explicit Matrix(Dim size, const Point offset, const Matrix<T> &parent);
    explicit Matrix(Dim size, const Matrix<T> &source);

    static Matrix<T> mul_classic(const Matrix<T> &lhs, const Matrix<T> &rhs);
    static Matrix<T> mul_winograd(const Matrix<T> &lhs, const Matrix<T> &rhs);
    static Matrix<T> mul_strassen(const Matrix<T> &lhs, const Matrix<T> &rhs);

    static void get_common_square_pow2(const Matrix<T> &src_a, const Matrix<T> &src_b, Matrix<T> &res_a, Matrix<T> &res_b);

    size_t get_index_or_throw(const size_t index) const;
    size_t get_index_or_throw(const size_t row, const size_t col) const;

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
    Matrix<T> &operator=(Matrix<T> &&source) noexcept;
    Matrix<T> operator+(const Matrix<T> &other) const;
    Matrix<T> operator-(const Matrix<T> &other) const;
    Matrix<T> operator*(const Matrix<T> &other) const;
    Matrix<T> operator*(const T &value) const;
    T &operator[](const size_t index);
    const T &operator[](const size_t index) const;

    Matrix<T> submatrix(Dim size, Point offset);
    Matrix<T> get_square_pow2() const;

    bool is_square() const;
    bool is_pow2() const;
    Dim get_size() const;

    T &at(const Point index);
    const T &at(const Point index) const;
};

#pragma region constructors

template<typename T>
Matrix<T>::Matrix(const Dim size, const Point offset, const Matrix<T> &parent)
        : _size(size)
        , _offset(offset)
        , _offset_plain(_offset.row * _total_cols + _offset.col)
        , _total_cols(parent._total_cols)
        , _data(parent._data) { }

template<typename T>
Matrix<T>::Matrix(Dim size, const Matrix<T> &source) : Matrix(size) {
    Dim src_size = source.get_size();
    if (src_size.get_rows() > _size.get_rows() || src_size.get_cols() > _size.get_cols()) {
        throw std::invalid_argument("Inner matrix must be smaller!");
    }

    for(int i = 0; i < src_size.get_rows(); ++i) {
        for(int j = 0; j < src_size.get_cols(); ++j) {
            at(i, j) = source.at(i, j);
        }
    }
}

template<typename T>
Matrix<T>::Matrix(const Dim size)
        : _size(size)
        , _offset()
        , _offset_plain()
        , _total_cols(size.get_cols())
        , _data(std::allocator<T>().allocate(size.get_card())) { }

template<typename T>
Matrix<T>::Matrix(const Dim size, const T &default_value)
        : Matrix(size) {
    std::fill(_data.get(), _data.get() + size.get_card(), default_value);
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T> &other)
        : Matrix(other._size) {
    int card = _size.get_card();
    for(int i = 0; i < card; ++i) {
        (*this)[i] = other[i];
    }
}

#pragma endregion

#pragma region operators

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &source) {
    if (this != &source) {
        _size = source._size;
        _offset = Point();
        _offset_plain = 0;
        _total_cols = _size.get_cols();

        size_t card = _size.get_card();
        _data = std::make_unique<T[]>(card);
        for(int i = 0; i < card; ++i) {
            (*this)[i] = source[i];
        }
    }

    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&source) noexcept {
    if (this != &source) {
        _data.reset();

        _size = source._size;
        _offset = source._offset;
        _offset_plain = source._offset_plain;
        _total_cols = _size.get_cols();
        _data = std::move(source._data);

        source._size = Dim(0, 0);
        source._offset = Point();
        source._offset_plain = 0;
        source._total_cols = 0;
        source._data.reset();
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
    // to be implemented
    return Matrix{other.get_size()};
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T &value) const {
    return mul(*this, value);
}

template<typename T>
T &Matrix<T>::operator[](const size_t index) {
     return _data.get()[get_index_or_throw(index)];
}

template<typename T>
const T &Matrix<T>::operator[](const size_t index) const {
    return _data.get()[get_index_or_throw(index)];
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
            return mul_winograd(lhs, rhs);
        case MulType::STRASSEN:
            return mul_strassen(lhs, rhs);
    }

    return mul_classic(lhs, rhs);
}

template<typename T>
Matrix<T> Matrix<T>::mul(const Matrix<T> &obj, const T &value) {
    Dim size = obj.get_size();
    Matrix<T> result{size};
    for(int i = 0; i < size.get_rows(); ++i) {
        for(int j = 0; j < size.get_cols(); ++j) {
            result[i][j] = obj[i][j] * value;
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::mul_classic(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    check_multiplicity(lhs, rhs);

    const int &m = lhs.get_size().get_rows();
    const int &n = lhs.get_size().get_cols();
    const int &k = rhs.get_size().get_cols();
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

    const size_t &m = lhs.get_size().get_rows();
    const size_t &n = lhs.get_size().get_cols();
    const size_t &k = rhs.get_size().get_cols();
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
    Matrix<T> &lhs_sqr = lhs;
    Matrix<T> &rhs_sqr = rhs;

    if (!lhs_sqr.is_pow2() || !rhs_sqr.is_pow2() || lhs.get_size() != lhs.get_size()) {
        get_common_square_pow2(lhs, rhs, lhs_sqr, rhs_sqr);
    }

    auto size = lhs_sqr.get_size();

    if (size == 2) {
        return mul_classic(lhs_sqr, rhs_sqr);
    }

    auto half = size >> 1;

    Matrix<T> res {size};
    Matrix<T> lhs_subs[4], rhs_subs[4], res_subs[4];
    for(int i = 0; i < size; i += half) {
        for(int j = 0; j < size; j += half) {
            lhs_subs[i] = lhs_sqr.submatrix(half, Point(i, j));
            rhs_subs[i] = rhs_sqr.submatrix(half, Point(i, j));
            res_subs[i] = res.submatrix(half, Point(i, j));
        }
    }

    auto p1 = mul_strassen((rhs_subs[1] - rhs_subs[3]), lhs_subs[0]);
    auto p2 = mul_strassen((lhs_subs[0] + lhs_subs[1]), rhs_subs[3]);
    auto p3 = mul_strassen((lhs_subs[1] + lhs_subs[2]), rhs_subs[0]);
    auto p4 = mul_strassen((rhs_subs[2] - rhs_subs[0]), lhs_subs[3]);
    auto p5 = mul_strassen((lhs_subs[0] + lhs_subs[3]), (rhs_subs[0] + rhs_subs[3]));
    auto p6 = mul_strassen((lhs_subs[3] - lhs_subs[1]), (rhs_subs[2] + rhs_subs[3]));
    auto p7 = mul_strassen((lhs_subs[0] - lhs_subs[2]), (rhs_subs[0] + rhs_subs[1]));

    res[0] = p4 + p5 + p6 - p2;
    res[1] = p1 + p2;
    res[2] = p3 + p4;
    res[3] = p1 + p5 - p3 - p7;

    return res;
}

#pragma endregion

#pragma region arithmetic

template<typename T>
void Matrix<T>::add(const Matrix<T> &rhs) {
    check_same_size(*this, rhs);

    size_t card = _size.get_card();
    for(int i = 0; i < card; ++i) {
        (*this)[i] += rhs[i];
    }
}

template<typename T>
void Matrix<T>::ddt(const Matrix<T> &rhs) {
    check_same_size(*this, rhs);

    size_t card = _size.get_card();
    for(int i = 0; i < card; ++i) {
        (*this)[i] -= rhs[i];
    }
}

template<typename T>
void Matrix<T>::mul(const T &value) {
    size_t card = _size.get_card();
    for(int i = 0; i < card; ++i) {
        (this)[i] *= value;
    }
}

#pragma endregion

#pragma region other

template<typename T>
void Matrix<T>::check_same_size(const Matrix<T> &a, const Matrix<T> &b) {
    if (a.get_size() != a.get_size()) {
        throw std::invalid_argument("Matrices must have the same dimensions!");
    }
}

template<typename T>
void Matrix<T>::check_multiplicity(const Matrix<T> &lhs, const Matrix<T> &rhs) {
    if (lhs.get_size().get_cols() != rhs.get_size().get_rows()) {
        throw std::invalid_argument("Dimensions don't match for multiplication!");
    }
}

template<typename T>
size_t Matrix<T>::get_index_or_throw(const size_t index) const {
    if (index > _size.get_card()) {
        throw std::out_of_range("Index out of Matrix bounds!");
    }

    if (_size.get_cols() != _total_cols) {
        size_t row = index / _total_cols;
        size_t col = index % _total_cols;
        return _offset_plain + col + row * _size.get_cols();
    }
    return index;
}

template<typename T>
size_t Matrix<T>::get_index_or_throw(const size_t row, const size_t col) const {
    if (row < _size.get_rows() && col < _size.get_cols()) {
        throw std::out_of_range("Index out of Matrix bounds!");
    }
    return _offset_plain + col + row * _size.get_cols();
}

template<typename T>
void
Matrix<T>::get_common_square_pow2(const Matrix<T> &src_a, const Matrix<T> &src_b, Matrix<T> &res_a, Matrix<T> &res_b) {
    if (src_a.is_pow2() && src_b.is_pow2() && src_a.get_size() == src_b.get_size()) {
        res_a = src_a;
        res_b = src_b;
        return;
    }
    size_t max_dim_a = std::max(src_a.get_size().get_rows(), src_a.get_size().get_cols());
    size_t max_dim_b = std::max(src_b.get_size().get_rows(), src_b.get_size().get_cols());
    size_t max_dim = std::max(max_dim_a, max_dim_b);
    size_t new_dim = Utils::bit_ceil(max_dim);
    res_a = Matrix<T>(Dim(new_dim, new_dim), src_a);
    res_b = Matrix<T>(Dim(new_dim, new_dim), src_b);
}

template<typename T>
Matrix<T> Matrix<T>::submatrix(Dim size, Point offset) {
    if (offset.row + size.get_rows() >= _size.get_rows() || offset.col + size.get_cols() > _size.get_cols()) {
        throw std::invalid_argument("Submatrix exceeds parent matrix dimensions!");
    }
    return Matrix(size, offset, *this);
}

template<typename T>
Matrix<T> Matrix<T>::get_square_pow2() const {
    if (is_pow2()) {
        return *this;
    }

    size_t max_dim = std::max(_size.get_rows(), _size.get_cols());
    size_t new_dim = Utils::bit_ceil(max_dim);
    return Matrix {Dim(new_dim, new_dim), *this};
}

template<typename T>
bool Matrix<T>::is_square() const {
    return _size.get_rows() == _size.get_cols();
}

template<typename T>
bool Matrix<T>::is_pow2() const {
    return _size.get_rows() == _size.get_cols() && Utils::is_power_of_two(_size.get_rows());
}

template<typename T>
Dim Matrix<T>::get_size() const {
    return _size;
}

template<typename T>
T &Matrix<T>::at(const Point index) {
    return _data.get()[get_index_or_throw(index.row, index.col)];
}

template<typename T>
const T &Matrix<T>::at(const Point index) const {
    return _data.get()[get_index_or_throw(index.row, index.col)];
}

#pragma endregion

#endif //LINAL_MATRIX_HPP