#include "linal/point.hpp"

Point::Point() : row(0), col(0) { }

Point::Point(const size_t row, const size_t col)
    : row(row)
    , col(col) { }

Point &Point::operator=(const Point &source) {
    if (this != &source) {
        row = source.row;
        col = source.col;
    }

    return *this;
}

