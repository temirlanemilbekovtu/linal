#pragma once

#ifndef LINAL_POINT_HPP
#define LINAL_POINT_HPP

#include "export.hpp"
#include "cstddef"

struct LINAL_API Point {
    size_t row;
    size_t col;

    explicit Point();
    explicit Point(const size_t row, const size_t col);

    Point &operator=(const Point &source);
};

#endif //LINAL_POINT_HPP
