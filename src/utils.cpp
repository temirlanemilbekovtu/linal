#include "linal/utils.hpp"

bool Utils::is_power_of_two(const uint32_t n)  {
    return (n & (n - 1)) == 0;
}

uint32_t Utils::bit_ceil(uint32_t n) {
    if (n <= 1) {
        return 1;
    }

    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    return n + 1;
}
