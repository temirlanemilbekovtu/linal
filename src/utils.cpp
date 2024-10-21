#include "linal/utils.hpp"

bool Utils::is_power_of_two(const int n)  {
    return n > 0 && (n & (n - 1)) == 0;
}

int Utils::bit_ceil(const int n) {
    int i = sizeof(n) * 8;
    bool bit;
    do {
        bit = (n >> --i) & 1;
    } while (bit != 1);
    return i + 2;
}
