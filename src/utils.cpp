#include "linal/utils.hpp"

bool Utils::is_power_of_two(const uint32_t n)  {
    return n > 0 && (n & (n - 1)) == 0;
}

uint32_t Utils::bit_ceil(uint32_t n) {
    int i = sizeof(n) * 8;
    bool bit;
    do {
        bit = (n >> --i) & 1;
    } while (bit != 1);
    return i + 2;
}
