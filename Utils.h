//==============================================================================
// Created on 2024/3/1.
// Copyright (c) WiwilZ. All rights reserved.
//==============================================================================

#pragma once


#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <bit>
#include <limits>


struct int128 {
    int64_t low;
    int64_t high;
};

struct uint128 {
    uint64_t low;
    uint64_t high;
};

template <typename T>
struct AddResult {
    T sum;
    bool carry;
};

template <typename T>
struct SubResult {
    T diff;
    bool borrow;
};

template <typename T>
struct DivResult {
    T quotient;
    T remainder;
};


namespace Naive {
    constexpr AddResult<uint64_t> Add(uint64_t a, uint64_t b, bool carry) noexcept {
        const uint64_t a0 = a + carry;
        carry = a > std::numeric_limits<uint64_t>::max() - carry;
        return {a0 + b, carry || a0 > std::numeric_limits<uint64_t>::max() - b};
    }

    constexpr SubResult<uint64_t> Sub(uint64_t a, uint64_t b, bool borrow) noexcept {
        const uint64_t a0 = a - borrow;
        borrow = a < static_cast<uint64_t>(borrow);
        return {a0 - b, borrow || a0 < b};
    }

    constexpr uint64_t MulHigh(uint64_t a, uint64_t b) noexcept {
        constexpr uint64_t mask = 0xffffffff;

        const uint64_t a1 = a >> 32;
        const uint64_t a0 = a & mask;
        const uint64_t b1 = b >> 32;
        const uint64_t b0 = b & mask;
        /*
         *   |    |     |  a0 * b0  |
         * + |    |  a1 * b0  |     |
         * - |    |  a0 * b1  |     |
         * + |  a1 * b1  |
         * --------------------------
         *   |   high   |          |
         */
        uint64_t t = (a0 * b0) >> 32;

        t += a1 * b0;
        uint64_t high = t >> 32;
        t &= mask;

        t += a0 * b1;
        high += t >> 32;

        high += a1 * b1;

        return high;
    }

    constexpr uint128 Mul(uint64_t a, uint64_t b) noexcept {
        constexpr uint64_t mask = 0xffffffff;

        const uint64_t a1 = a >> 32;
        const uint64_t a0 = a & mask;
        const uint64_t b1 = b >> 32;
        const uint64_t b0 = b & mask;
        /*
         *   |    |     |  a0 * b0  |
         * + |    |  a1 * b0  |     |
         * - |    |  a0 * b1  |     |
         * + |  a1 * b1  |
         * --------------------------
         *   |   high   |   low     |
         */
        uint64_t low = a0 * b0;
        uint64_t t = low >> 32;
        low &= mask;

        t += a1 * b0;
        uint64_t high = t >> 32;
        t &= mask;

        t += a0 * b1;
        low |= t << 32;
        high += t >> 32;

        high += a1 * b1;

        return {low, high};
    }

    constexpr DivResult<uint64_t> Div(uint128 a, uint64_t b) noexcept {
        constexpr uint64_t mask = 0xffffffff;

        uint64_t ah = a.high;
        uint64_t al = a.low;

        if (ah == 0) {
            return {ah / b, al % b};
        }

        const unsigned shift = std::countl_zero(b);
        if (shift > 0) {
            b <<= shift;
            ah = (ah << shift) | (al >> (64 - shift));
            al <<= shift;
        }

        /*
         * b的MSB为1 => qhat - 2 <= q <= qhat
         * [a2 a1 a0]
         *    [b1 b0]
         * q = [a21 a0] / b
         * estimate qhat = a21 / b1 ... r  =>  a21 = q * b1 + r
         * qhat * b = q * b1 * 2^32 + q * b0
         *        a = a21 * 2^32 + a0
         *          = q * b1 * 2^32 + r * 2^32 + a0
         * t0 = qhat * b0
         * t1 = r * 2^32 + a0
         * t0 - t1 >  b => qhat * b - a >  b => (qhat - 1) * b >  a => q = qhat - 2
         * t0 - t1 >  0 => qhat * b - a >  0 =>       qhat * b >  a => q = qhat - 1
         * t0 - t1 <= 0 => qhat * b - a <= 0 =>       qhat * b <= a => q = qhat
        */

        const uint64_t al1 = al >> 32;
        const uint64_t al0 = al & mask;
        const uint64_t b1 = b >> 32;
        const uint64_t b0 = b & mask;

        uint64_t q1 = ah / b1;
        uint64_t r = ah % b1;
        uint64_t t0 = q1 * b0;
        uint64_t t1 = (r << 32) | al1;
        if (t0 > t1) {
            q1 -= ((t0 - t1) > b) + 1;
        }
        uint64_t u = ((ah << 32) | al1) - q1 * b; // remainder is guaranteed to be 64 bits

        uint64_t q0 = u / b1;
        r = u % b1;
        t0 = q0 * b0;
        t1 = (r << 32) | al0;
        if (t0 > t1) {
            q0 -= ((t0 - t1) > b) + 1;
        }
        u = ((u << 32) | al0) - q0 * b;

        return {(q1 << 32) | q0, u >> shift};
    }
}


namespace Utils {
    constexpr int GetSign(std::integral auto v) noexcept {
        return (v > 0) - (v < 0);
    }

    template <std::integral T>
    constexpr T Abs(T v) noexcept {
        return v < 0 ? -static_cast<std::make_unsigned_t<T>>(v) : v;
    }

    constexpr uint64_t* Extend(uint64_t* const old_ptr, size_t old_size, size_t new_size) {
        if (new_size <= old_size) {
            return old_ptr;
        }
        uint64_t* ptr{};
        if !consteval {
            ptr = static_cast<uint64_t*>(std::realloc(old_ptr, new_size * sizeof(uint64_t)));
        }
        if (ptr == nullptr) {
            ptr = new uint64_t[new_size];
            std::copy_n(old_ptr, old_size, ptr);
            delete[] old_ptr;
        }
        return ptr;
    }

    constexpr std::strong_ordering Compare(const uint64_t* const a, const uint64_t* const b, int64_t size) noexcept {
        for (int64_t i = size - 1; i >= 0; --i) {
            if (const auto cmp = a[i] <=> b[i]; cmp != 0) {
                return cmp;
            }
        }
        return std::strong_ordering::equal;
    }

    constexpr int64_t Normalize(const uint64_t* const data, int64_t size) noexcept {
        int64_t i = size;
        while (i > 0 && data[i - 1] == 0) {
            --i;
        }
        return i;
    }

    constexpr bool Increase(uint64_t* const data, int64_t size) noexcept {
        // 最低位开始，连续的uint64_t::max结果为0，向第一个非uint64_t::max位x进1
        // x += 1，不向前进位
        int64_t i = 0;
        while (i < size && data[i] == std::numeric_limits<uint64_t>::max()) {
            data[i] = 0;
            ++i;
        }
        if (i < size) {
            ++data[i];
            return false;
        }
        return true;
    }

    constexpr bool Increase(uint64_t* const dst, const uint64_t* const src, int64_t size) noexcept {
        int64_t i = 0;
        while (i < size && src[i] == std::numeric_limits<uint64_t>::max()) {
            dst[i] = 0;
            ++i;
        }
        if (i < size) {
            dst[i] = src[i] + 1;
            std::copy(src + i + 1, src + size, dst + i + 1);
            return false;
        }
        return true;
    }

    constexpr void Decrease(uint64_t* const data, int64_t size) noexcept {
        // data最高位不为0
        // 最低位开始，连续的0结果为uint64_t::max，向第一个非0位x借1
        // x -= 1，不向前借位
        int64_t i = 0;
        while (data[i] == 0) {
            data[i] = std::numeric_limits<uint64_t>::max();
            ++i;
        }
        --data[i];
    }

    constexpr void Decrease(uint64_t* const dst, const uint64_t* const src, int64_t size) noexcept {
        int64_t i = 0;
        while (src[i] == 0) {
            dst[i] = std::numeric_limits<uint64_t>::max();
            ++i;
        }
        dst[i] = src[i] - 1;
        std::copy(src + i + 1, src + size, dst + i + 1);
    }


    constexpr AddResult<uint64_t> Add(uint64_t a, uint64_t b, bool carry) noexcept {
#if defined(__SIZEOF_INT128__)
        const auto result = static_cast<__uint128_t>(a) + b + carry;
        return {static_cast<uint64_t>(result), static_cast<bool>(result >> 64)};
#elif defined(_MSC_VER) && defined(_M_X64) && !defined(_M_ARM64EC)
        if consteval {
            return Naive::Add(a, b, carry);
        } else {
            uint64_t result;
            carry = _addcarryx_u64(carry, a, b, &result);
            return {result, carry};
        }
#else
        return Naive::Add(a, b, carry);
#endif
    }

    constexpr SubResult<uint64_t> Sub(uint64_t a, uint64_t b, bool borrow) noexcept {
#if defined(__SIZEOF_INT128__)
        const auto result = static_cast<__uint128_t>(a) - b - borrow;
        return {static_cast<uint64_t>(result), static_cast<bool>(result >> 64)};
#elif defined(_MSC_VER) && defined(_M_X64) && !defined(_M_ARM64EC)
        if consteval {
            return Naive::Sub(a, b, carry);
        } else {
            uint64_t result;
            borrow = _subborrow_u64(borrow, a, b, &result);
            return {result, borrow};
        }
#else
        return Naive::Sub(a, b, carry);
#endif
    }

    constexpr uint64_t MulHigh(uint64_t a, uint64_t b) noexcept {
#if defined(__SIZEOF_INT128__)
        return (static_cast<__uint128_t>(a) * b) >> 64;
#elif defined(_MSC_VER) && defined(_M_X64) && !defined(_M_ARM64EC)
        if consteval {
            return Naive::MulHigh(a, b);
        } else {
            return __umulh(a, b);
        }
#else
        return Naive::MulHigh(a, b);
#endif
    }

    constexpr uint128 Mul(uint64_t a, uint64_t b) noexcept {
#if defined(__SIZEOF_INT128__)
        /*
        // mulx 积的高位(dst1), 积的低位(dst2), 乘数(src)
        // 另一乘数为rdx
        uint128 result;
        __asm__("mulx %[a], %%rax, %%rdx"
                : "=d"(result.high), "=a"(result.low)
                : [a] "r"(a), "d"(b));
        return result;
        */
        const __uint128_t result = static_cast<__uint128_t>(a) * b;
        return {static_cast<uint64_t>(result), static_cast<uint64_t>(result >> 64)};
#elif defined(_MSC_VER) && defined(_M_X64) && !defined(_M_ARM64EC)
        if consteval {
            return Naive::Mul(a, b);
        } else {
            uint64_t high;
            // const uint64_t low = _umul128(a, b, &high);
            const uint64_t low = _mulx_u64(a, b, &high);
            return {low, high};
        }
#else
        return Naive::Mul(a, b);
#endif
    }

    constexpr DivResult<uint64_t> Div(uint128 a, uint64_t b) noexcept {
        assert(a.high < b && b != 0);

#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__)
        if consteval {
            return Naive::Div(a, b);
        } else {
            DivResult<uint64_t> result;
            __asm__("divq %[b]"
                : "=d"(result.remainder), "=a"(result.quotient)
                : "d"(a.high), "a"(a.low), [b] "r"(b));
            return result;
        }
#elif defined(_MSC_VER) && defined(_M_X64) && !defined(_M_ARM64EC)
        if consteval {
            return Naive::Div(a, b);
        } else {
            uint64_t remainder;
            const uint64_t quotient = _udiv128(a.high, a.low, b, &remainder);
            return {quotient, remainder};
        }
#else
        return Naive::Div(a, b);
#endif
    }
}




// #if defined(_MSC_VER) && !defined(__clang__)
// // 寄存器传参顺序：rcx rdx r8 r9
//
// int64_t MulHigh(int64_t a, int64_t b) {
//     return __mulh(a, b);
// }
//
// int128 Mul(int64_t a, int64_t b) {
//     int64_t high;
//     const int64_t low = _mul128(a, b, &high);
//     return {low, high};
// }
//
// DivResult<int64_t> Div(int64_t ah, int64_t al, int64_t b) {
//     int64_t remainder;
//     const int64_t quotient = _div128(ah, al, b, &remainder);
//     return {quotient, remainder};
// }
//
// #else
// // 寄存器传参顺序：rdi rsi rdx rcx r8 r9
//
// int64_t MulHigh(int64_t a, int64_t b) {
//     return (static_cast<__int128_t>(a) * b) >> 64;
// }
//
// int128 Mul(int64_t a, int64_t b) {
//     /*
//     // imul 乘数
//     // 另一乘数为rax
//     // 结果为rdx:rax
//     int128 result;
//     __asm__("imulq %[a]"
//             : "=d"(result.high), "=a"(result.low)
//             : [a] "r"(a), "a"(b)
//             : "cc");
//     return result;
//     */
//     const auto result = static_cast<__int128_t>(a) * b;
//     return {result, result >> 64};
// }
//
// DivResult<int64_t> Div(int64_t ah, int64_t al, int64_t b) {
//     // idiv 除数
//     // 被除数为rdx:rax
//     // 余数为rdx，商为rax
//     DivResult<int64_t> result;
//     __asm__("idivq %[b]"
//             : "=d"(result.remainder), "=a"(result.quotient)
//             : "d"(ah), "a"(al), [b] "r"(b));
//     return result;
// }
// #endif
