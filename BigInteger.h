//==============================================================================
// Created on 24-3-2.
// Copyright (c) WiwilZ. All rights reserved.
//==============================================================================

#pragma once

#include "Utils.h"
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <bit>
#include <compare>
#include <concepts>
#include <iostream>
#include <limits>
#include <ranges>
#include <span>
#include <string>
#include <string_view>


class BigInteger {
private:
    // number = sign(size) * data
    int64_t size{};
    uint64_t* data{};

    static constexpr uint64_t MaxDigit = std::numeric_limits<uint64_t>::max();

    constexpr int64_t GetLength() const noexcept {
        const auto x = Utils::Abs(size);
        return x | (x == 0);
    }

    constexpr BigInteger(int64_t size, uint64_t* data) noexcept : size(size), data(data) {}

    template <std::unsigned_integral T>
    static constexpr BigInteger FromNegative(T v) {
        if (v > MaxDigit) {
            return BigInteger{-2, new uint64_t[2]{static_cast<uint64_t>(v), static_cast<uint64_t>(v >> 64)}};
        }
        return BigInteger{-(v != 0), new uint64_t[1]{static_cast<uint64_t>(v)}};
    }

public:
    static const BigInteger zero;
    static const BigInteger one;
    static const BigInteger minus_one;

public:// Big 5
    constexpr BigInteger() noexcept = default;

    constexpr BigInteger(const BigInteger& v) : size(v.size) {
        const auto len = Utils::Abs(v.size);
        data = new uint64_t[len];
        std::copy_n(v.data, len, data);
    }

    constexpr BigInteger& operator=(const BigInteger& v) noexcept {
        if (this != &v) {
            const auto len = Utils::Abs(v.size);
            data = Utils::Extend(data, Utils::Abs(size), len);
            std::copy_n(v.data, len, data);
            size = v.size;
        }
        return *this;
    }

    constexpr BigInteger(BigInteger&& v) noexcept {
        size = std::exchange(v.size, 0);
        data = std::exchange(v.data, nullptr);
    }

    constexpr BigInteger& operator=(BigInteger&& v) noexcept {
        if (this != &v) {
            delete[] data;
            size = std::exchange(v.size, 0);
            data = std::exchange(v.data, nullptr);
        }
        return *this;
    }

    constexpr ~BigInteger() noexcept {
        delete[] data;
    }

public:
    template <std::integral T>
    constexpr explicit BigInteger(T v) {
        const std::make_unsigned_t<T> u = Utils::Abs(v);
        if (u > std::numeric_limits<uint64_t>::max()) {
            size = 2;
            data = new uint64_t[2];
            data[1] = u >> 64;
        } else {
            size = u != 0;
            data = new uint64_t[1];
        }
        data[0] = u;
        if (v < 0) {
            size = -size;
        }
    }

    template <std::integral T>
    constexpr BigInteger& operator=(T v) {
        const std::make_unsigned_t<T> u = Utils::Abs(v);
        if (u > std::numeric_limits<uint64_t>::max()) {
            data = Utils::Extend(data, GetLength(), 2);
            size = 2;
            data[1] = u >> 64;
        } else {
            size = u != 0;
        }
        data[0] = u;
        if (v < 0) {
            size = -size;
        }
        return *this;
    }

    constexpr explicit BigInteger(const std::string& str, int radix = 10) {
        if (radix < 2 || radix > 36) {
            throw std::out_of_range("radix must be between 2 and 36");
        }
        //TODO
    }

    constexpr BigInteger& operator=(const std::string& val) {
        //TODO
        return *this;
    }

public:
    constexpr int sign() const noexcept {
        return Utils::GetSign(size);
    }

    constexpr bool is_zero() const noexcept {
        return size == 0;
    }

    constexpr bool is_positive() const noexcept {
        return size > 0;
    }

    constexpr bool is_negative() const noexcept {
        return size < 0;
    }

    constexpr bool is_even() const noexcept {
        return (data[0] & 1) == 0;
    }

    constexpr bool is_odd() const noexcept {
        return (data[0] & 1) == 1;
    }

    constexpr size_t bit_count() const noexcept {
        if (size < 0) {
            return -1;
        }
        if (size == 0) {
            return 0;
        }
        const size_t len = size - 1;
        return len * 64 + std::bit_width(data[len]);
    }

    constexpr size_t leading_zeros() const noexcept {
        if (size < 0) {
            return 0;
        }
        if (size == 0) {
            return -1;
        }
        return std::countl_zero(data[size - 1]);
    }

    constexpr size_t tailing_zeros() const noexcept {
        if (size == 0) {
            return -1;
        }
        /*
         * -a = ~a + 1
         * 最低位开始，连续的0结果仍为0，向第一个非0位x进1
         * ~x + 1 = -x 不为0
         */
        size_t i = 0;
        while (data[i] == 0) {
            ++i;
        }
        const uint64_t x = size < 0 ? -data[i] : data[i];
        return i * 64 + std::countr_zero(x);
    }

    constexpr bool is_power_of_2() const noexcept {
        if (size <= 0) {
            return false;
        }
        for (auto* p = data; p != data + size - 1; ++p) {
            if (*p != 0) {
                return false;
            }
        }
        return std::has_single_bit(data[size - 1]);
    }

    std::string to_string(int radix = 10) const {
        if (radix < 2 || radix > 36) {
            throw std::out_of_range("radix must be between 2 and 36");
        }
        if (is_zero()) {
            return "0";
        }
        //TODO
        return "";
    }

public:// 强制类型转换运算符
    constexpr explicit operator bool() const noexcept {
        return !is_zero();
    }

    template <std::integral T>
    constexpr explicit operator T() const noexcept {
        const std::make_unsigned_t<T> u = data[0];
        return static_cast<T>(size < 0 ? -u : u);
    }

#ifdef __SIZEOF_INT128__
    constexpr explicit operator __uint128_t() const noexcept {
        const __uint128_t h = size > 1 || size < -1 ? data[1] : 0;
        const __uint128_t u = (h << 64) | data[0];
        return size < 0 ? -u : u;
    }

    constexpr explicit operator __int128_t() const noexcept {
        return static_cast<__int128_t>(static_cast<__uint128_t>(*this));
    }
#endif

    constexpr explicit operator double() const noexcept {
        //TODO
        return 0.0;
    }

    constexpr explicit operator float() const noexcept {
        return static_cast<float>(static_cast<double>(*this));
    }

public:// 比较运算
    friend constexpr std::strong_ordering operator<=>(const BigInteger& lhs, const BigInteger& rhs) noexcept {
        auto cmp = lhs.size <=> rhs.size;
        if (cmp != 0) {
            return cmp;
        }
        cmp = Utils::Compare(lhs.data, rhs.data, lhs.GetLength());
        return lhs.size < 0 ? 0 <=> cmp : cmp;
    }

    friend constexpr bool operator==(const BigInteger& lhs, const BigInteger& rhs) noexcept {
        return lhs.size == rhs.size && Utils::Compare(lhs.data, rhs.data, lhs.GetLength()) == 0;
    }

    friend constexpr std::strong_ordering operator<=>(const BigInteger& lhs, std::integral auto rhs) noexcept {
        auto cmp = lhs.size <=> Utils::GetSign(rhs);
        if (cmp != 0) {
            return cmp;
        }
        const uint64_t u = lhs.data[0];
        const uint64_t v = Utils::Abs(rhs);
        return rhs < 0 ? v <=> u : u <=> v;
    }

    friend constexpr bool operator==(const BigInteger& lhs, std::integral auto rhs) noexcept {
        return lhs.size == Utils::GetSign(rhs) && lhs.data[0] == static_cast<uint64_t>(Utils::Abs(rhs));
    }

#ifdef __SIZEOF_INT128__
    friend constexpr std::strong_ordering operator<=>(const BigInteger& lhs, __uint128_t rhs) noexcept {
        if (lhs.size < 0) {
            return std::strong_ordering::less;
        }
        if (lhs.size > 2) {
            return std::strong_ordering::greater;
        }
        const __uint128_t h = lhs.size == 2 ? lhs.data[1] : 0;
        const __uint128_t u = (h << 64) | lhs.data[0];
        return u <=> rhs;
    }

    friend constexpr bool operator==(const BigInteger& lhs, __uint128_t rhs) noexcept {
        if (lhs.size < 0 || lhs.size > 2) {
            return false;
        }
        const __uint128_t h = lhs.size == 2 ? lhs.data[1] : 0;
        const __uint128_t u = (h << 64) | lhs.data[0];
        return u == rhs;
    }

    friend constexpr std::strong_ordering operator<=>(const BigInteger& lhs, __int128_t rhs) noexcept {
        if (lhs.size < -2) {
            return std::strong_ordering::less;
        }
        if (lhs.size > 2) {
            return std::strong_ordering::greater;
        }
        const __uint128_t h = lhs.size == -2 || lhs.size == 2 ? lhs.data[1] : 0;
        const __uint128_t u = (h << 64) | lhs.data[0];
        const __uint128_t v = Utils::Abs(rhs);
        return rhs < 0 ? v <=> u : u <=> v;
    }

    friend constexpr bool operator==(const BigInteger& lhs, __int128_t rhs) noexcept {
        if (lhs.size < -2 || lhs.size > 2) {
            return false;
        }
        const __uint128_t h = lhs.size == -2 || lhs.size == 2 ? lhs.data[1] : 0;
        const __uint128_t u = (h << 64) | lhs.data[0];
        const __uint128_t v = Utils::Abs(rhs);
        return u == v;
    }
#endif

    friend constexpr std::strong_ordering operator<=>(std::integral auto lhs, const BigInteger& rhs) noexcept {
        return 0 <=> (rhs <=> lhs);
    }

    friend constexpr bool operator==(std::integral auto lhs, const BigInteger& rhs) noexcept {
        return rhs == lhs;
    }

public:// 单目运算
    constexpr BigInteger operator~() const {
        //  x: ~<x> = -x - 1 = -<x + 1>
        // -x: ~<-x> = ~[~(x - 1)] = <x - 1>
        int64_t res_size;
        auto* const res_data = new uint64_t[GetLength() + 1];
        if (size < 0) {
            res_size = -size;
            Utils::Decrease(res_data, data, res_size);
            res_size -= res_data[res_size - 1] == 0;
        } else {
            res_size = size;
            const auto carry = Utils::Increase(res_data, data, res_size);
            res_data[res_size] = carry;
            res_size += carry;
            res_size = -res_size;
        }
        return BigInteger{res_size, res_data};
    }

    constexpr void to_bit_not() {
        if (size < 0) {
            size = -size;
            Utils::Decrease(data, size);
            size -= data[size - 1] == 0;
        } else {
            if (const auto carry = Utils::Increase(data, size); carry) {
                data = Utils::Extend(data, size, size + 1);
                data[size] = 1;
                ++size;
            }
            size = -size;
        }
    }

    constexpr BigInteger operator+() const noexcept {
        return *this;
    }

    constexpr BigInteger operator-() const {
        const auto res_len = GetLength();
        auto* const res_data = new uint64_t[res_len];
        std::copy_n(data, res_len, res_data);
        return BigInteger{-size, res_data};
    }

    constexpr void to_negative() noexcept {
        size = -size;
    }

    constexpr BigInteger abs() const {
        const auto res_len = GetLength();
        auto* const res_data = new uint64_t[res_len];
        std::copy_n(data, res_len, res_data);
        return BigInteger{size == 0 ? 0 : res_len, res_data};
    }

    constexpr void to_abs() noexcept {
        if (size < 0) {
            size = -size;
        }
    }

    constexpr BigInteger& operator++() {
        //  x: ++<x> = <x + 1>
        // -x: ++<-x> = -<x> + 1 = -<x - 1>
        if (size < 0) {
            size = -size;
            Utils::Decrease(data, size);
            size -= data[size - 1] == 0;
            size = -size;
        } else if (const auto carry = Utils::Increase(data, size); carry) {
            data = Utils::Extend(data, size, size + 1);
            data[size] = 1;
            ++size;
        }
        return *this;
    }

    constexpr BigInteger operator++(int) {
        BigInteger res(*this);
        ++*this;
        return res;
    }

    constexpr BigInteger& operator--() {
        //  x: --<x> = <x - 1>
        // -x: --<-x> = <-x> - 1 = -<x + 1>
        if (size > 0) {
            Utils::Decrease(data, size);
            size -= data[size - 1] == 0;
        } else {
            size = -size;
            if (const auto carry = Utils::Increase(data, size); carry) {
                data = Utils::Extend(data, size, size + 1);
                data[size] = 1;
                ++size;
            }
            size = -size;
        }
        return *this;
    }

    constexpr BigInteger operator--(int) {
        BigInteger res(*this);
        --*this;
        return res;
    }

public:// & | ^ 运算
    friend constexpr BigInteger operator&(const BigInteger& lhs, const BigInteger& rhs) {
        if (lhs.size == 0 || rhs.size == 0) {
            return zero;
        }

        const BigInteger& u = lhs.size >= rhs.size ? lhs : rhs;
        const BigInteger& v = &u == &lhs ? rhs : lhs;

        // u >= v > 0  结果为正数或0，长度 <= v.size
        if (v.size > 0) {
            auto res_size = v.size;
            while (res_size > 0 && (u.data[res_size - 1] & v.data[res_size - 1]) == 0) {
                --res_size;
            }
            if (res_size == 0) {
                return zero;
            }
            auto* const res_data = new uint64_t[res_size];
            for (size_t i = 0; i < res_size; ++i) {
                res_data[i] = u.data[i] & v.data[i];
            }
            return BigInteger{res_size, res_data};
        }

        const auto vlen = -v.size;
        // 0 > u >= v  结果为负数，长度 <= vlen + 1
        if (u.size < 0) {
            /*
             *   -(<-u> & <-v>)
             * = -[~(u - 1) & ~(v - 1)]
             * = ~[~(u - 1) & ~(v - 1)] + 1
             * = [(u - 1) | (v - 1)] + 1
             */
            const auto ulen = -u.size;

            auto* const utmp = new uint64_t[ulen];
            Utils::Decrease(utmp, u.data, ulen);

            auto res_len = ulen;
            auto* const res_data = new uint64_t[res_len + 1];
            Utils::Decrease(res_data, v.data, vlen);

            for (size_t i = 0; i < ulen; ++i) {
                res_data[i] |= utmp[i];
            }

            delete[] utmp;
            res_len -= res_data[res_len - 1] == 0;

            const auto carry = Utils::Increase(res_data, res_len);
            res_data[res_len] = carry;
            res_len += carry;

            return BigInteger{-res_len, res_data};
        }

        // u > 0 > v  结果为正数，长度 <= u.size
        // u & -v = u & ~(v - 1)
        auto res_size = u.size;
        auto* const res_data = new uint64_t[res_size];
        if (u.size <= vlen) {
            /*
             * 00 00 u1 u0
             * 11 v2 v1 v0
             *
             * 00 u2 u1 u0
             * 11 v2 v1 v0
             */
            Utils::Decrease(res_data, v.data, u.size);
            while (res_size > 0 && (u.data[res_size - 1] & ~res_data[res_size - 1]) == 0) {
                --res_size;
            }
            res_data[0] = u.data[0] & ~res_data[0];
            for (size_t i = 1; i < res_size; ++i) {
                res_data[i] = u.data[i] & ~res_data[i];
            }
        } else {
            /*
            * u3 u2 u1 u0
            * 11 11 v1 v0
            */
            Utils::Decrease(res_data, v.data, vlen);
            for (size_t i = 0; i < vlen; ++i) {
                res_data[i] = u.data[i] & ~res_data[i];
            }
            std::copy(u.data + vlen, u.data + u.size, res_data + vlen);
        }
        return BigInteger{res_size, res_data};
    }

    constexpr BigInteger& operator&=(const BigInteger& v) {
        if (this != &v) {
            *this = *this & v;
        }
        return *this;
    }

    template <std::integral T>
    friend constexpr BigInteger operator&(const BigInteger& lhs, T rhs) {
        return lhs & BigInteger(rhs);
    }

    template <std::integral T>
    constexpr BigInteger& operator&=(T v) {
        *this = *this & v;
        return *this;
    }

    template <std::integral T>
    friend constexpr BigInteger operator&(T lhs, const BigInteger& rhs) {
        return rhs & lhs;
    }

    friend constexpr BigInteger operator|(const BigInteger& lhs, const BigInteger& rhs) {
        if (lhs.size == 0) {
            return rhs;
        }
        if (rhs.size == 0) {
            return lhs;
        }

        const BigInteger& u = lhs.size >= rhs.size ? lhs : rhs;
        const BigInteger& v = &u == &lhs ? rhs : lhs;

        // u >= v > 0  结果为正数，长度 <= u.size
        if (v.size > 0) {
            auto* const res_data = new uint64_t[u.size];
            for (size_t i = 0; i < v.size; ++i) {
                res_data[i] = u.data[i] | v.data[i];
            }
            std::copy(u.data + v.size, u.data + u.size, res_data + v.size);
            return BigInteger{u.size, res_data};
        }

        const auto vlen = -v.size;
        // 0 > u >= v  结果为负数，长度 <= ulen
        if (u.size < 0) {
            /*
             *   -(<-u> | <-v>)
             * = -[~(u - 1) | ~(v - 1)]
             * = ~[~(u - 1) | ~(v - 1)] + 1
             * = [(u - 1) & (v - 1)] + 1
            */
            const auto ulen = -u.size;

            auto res_len = ulen;
            auto* const res_data = new uint64_t[res_len];
            Utils::Decrease(res_data, u.data, res_len);

            auto* const vtmp = new uint64_t[res_len];
            Utils::Decrease(vtmp, v.data, res_len);

            while (res_len > 0 && (res_data[res_len - 1] & vtmp[res_len - 1]) == 0) {
                --res_len;
            }
            res_data[0] &= vtmp[0];
            for (size_t i = 1; i < res_len; ++i) {
                res_data[i] &= vtmp[i];
            }
            delete[] vtmp;

            const auto carry = Utils::Increase(res_data, res_len);
            res_data[res_len] = carry;
            res_len += carry;

            return BigInteger{-res_len, res_data};
        }

        // u > 0 > v  结果为负数，长度 <= vlen
        /*
         *   -(<u> | <-v>)
         * = -[u | ~(v - 1)]
         * = ~[u | ~(v - 1)] + 1
         * = [~u & (v - 1)] + 1
        */
        auto res_len = vlen;
        auto* const res_data = new uint64_t[res_len];
        Utils::Decrease(res_data, v.data, vlen);

        if (u.size < vlen) {
            /*
             * 11 11 u1 u0
             * 00 v2 v1 v0
             */
            for (size_t i = 0; i < u.size; ++i) {
                res_data[i] &= ~u.data[i];
            }
        } else {
            /*
            * 11 u2 u1 u0
            * 00 00 v1 v0
            *
            * 11 u2 u1 u0
            * 00 v2 v1 v0
            */
            while (res_len > 0 && (~u.data[res_len - 1] & res_data[res_len - 1]) == 0) {
                --res_len;
            }
            res_data[0] = u.data[0] & ~res_data[0];
            for (size_t i = 1; i < res_len; ++i) {
                res_data[i] &= ~u.data[i];
            }
        }

        const auto carry = Utils::Increase(res_data, res_len);
        res_data[res_len] = carry;
        res_len += carry;
        return BigInteger{-res_len, res_data};
    }

    constexpr BigInteger& operator|=(const BigInteger& v) {
        if (this != &v) {
            *this = *this | v;
        }
        return *this;
    }

    template <std::integral T>
    friend constexpr BigInteger operator|(const BigInteger& lhs, T rhs) {
        return lhs | BigInteger(rhs);
    }

    template <std::integral T>
    constexpr BigInteger& operator|=(T v) {
        *this = *this | v;
        return *this;
    }

    template <std::integral T>
    friend constexpr BigInteger operator|(T lhs, const BigInteger& rhs) {
        return rhs | lhs;
    }

public:
};


const BigInteger BigInteger::zero{};
const BigInteger BigInteger::one{1, new uint64_t[1]{1}};
const BigInteger BigInteger::minus_one{-1, new uint64_t[1]{1}};
