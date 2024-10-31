#pragma once

#ifndef LINAL_UTILS_HPP
#define LINAL_UTILS_HPP

#include "cstdint"

class Utils {
public:
    static uint32_t bit_ceil(const uint32_t n);
    static bool is_power_of_two(const uint32_t n);
};

#endif //LINAL_UTILS_HPP
