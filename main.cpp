#include <iostream>
#include "Utils.h"



int main() {
    constexpr uint64_t a = 1;
    constexpr uint64_t b = 1;
    constexpr bool carry = true;
    constexpr auto result = Utils::Add(a, b, carry);
    std::cout << result.sum << " " << result.carry << std::endl;
    return 0;
}
