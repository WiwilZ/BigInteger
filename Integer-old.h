/*
 * Created by WiwilZ on 2022/6/28.
 */

#pragma once

#include <limits>
#include <sstream>
#include <regex>
#include <cmath>
#include <x86intrin.h>
#include <cassert>
#include <span>


class Integer {
public:
	template <std::integral T>
	static constexpr int sign(T v) noexcept {
		return (v > 0) - (v < 0);
	}

	static constexpr int sign(__int128_t v) noexcept {
		return (v > 0) - (v < 0);
	}

	static constexpr int sign(__uint128_t v) noexcept {
		return v > 0;
	}

	template <std::integral T>
	static constexpr std::make_unsigned_t<T> abs(T v) noexcept {
		return v < 0 ? -v : v;
	}

	static constexpr __uint128_t abs(__int128_t v) noexcept {
		return v < 0 ? -v : v;
	}

	static constexpr __uint128_t abs(__uint128_t v) noexcept {
		return v;
	}

	template <std::unsigned_integral T>
	static constexpr size_t leading_zeros(T v) noexcept {
		return std::countl_zero(v);
	}

	static constexpr size_t leading_zeros(__uint128_t v) noexcept {
		const auto low = static_cast<uint64_t>(v);
		const auto high = static_cast<uint64_t>(v >> 64);
		return high == 0 ? 64 + leading_zeros(low) : leading_zeros(high);
	}

	template <std::unsigned_integral T>
	static constexpr size_t tailing_zeros(T v) noexcept {
		return std::countr_zero(v);
	}

	static constexpr size_t tailing_zeros(__uint128_t v) noexcept {
		const auto low = static_cast<uint64_t>(v);
		const auto high = static_cast<uint64_t>(v >> 64);
		return low == 0 ? 64 + tailing_zeros(high) : tailing_zeros(low);
	}

	template <std::unsigned_integral T>
	static constexpr size_t bit_count(T v) noexcept {
		return std::bit_width(v);
	}

	static constexpr size_t bit_count(__uint128_t v) noexcept {
		const auto low = static_cast<uint64_t>(v);
		const auto high = static_cast<uint64_t>(v >> 64);
		auto res = bit_count(low);
		if (high != 0) {
			res += bit_count(high);
		}
		return res;
	}

	template <std::unsigned_integral T>
	static constexpr bool is_power_of_2(T v) noexcept {
		return std::has_single_bit(v);
	}

	static constexpr bool is_power_of_2(__uint128_t v) noexcept {
		const auto low = static_cast<uint64_t>(v);
		const auto high = static_cast<uint64_t>(v >> 64);
		return high == 0 ? is_power_of_2(low) : low == 0 && is_power_of_2(high);
	}

private:
	static constexpr uint64_t Add(unsigned long long x, unsigned long long y, bool& carry) noexcept {
		unsigned long long res;
		carry = _addcarryx_u64(carry, x, y, &res);
		return res;
	}

	static constexpr uint64_t Add(uint64_t x, bool& carry) noexcept {
		const uint64_t res = x + carry;
		carry = !res;
		return res;
	}

	static constexpr uint64_t Sub(unsigned long long x, unsigned long long y, bool& borrow) noexcept {
		unsigned long long res;
		borrow = _subborrow_u64(borrow, x, y, &res);
		return res;
	}

	static constexpr uint64_t Sub(uint64_t x, bool& borrow) noexcept {
		const uint64_t res = x - borrow;
		borrow = !x && borrow;
		return res;
	}

	static constexpr uint64_t MulAdd(uint64_t x, uint64_t y, uint64_t& carry) noexcept {
		const auto val = static_cast<__uint128_t>(x) * y + carry;
		carry = static_cast<uint64_t>(val >> 64);
		return static_cast<uint64_t>(val);
	}

	static constexpr uint64_t MulAdd(uint64_t x, uint64_t y, uint64_t z, uint64_t& carry) noexcept {
		const auto val = static_cast<__uint128_t>(x) * y + z + carry;
		carry = static_cast<uint64_t>(val >> 64);
		return static_cast<uint64_t>(val);
	}

	static constexpr uint64_t Div(uint64_t& high, uint64_t low, uint64_t y) {
		const auto x = static_cast<__uint128_t>(high) << 64 | low;
		const auto res = x / y;
		high = x - res * y;
		return res;
	}

	static constexpr uint64_t Complement(uint64_t val, uint64_t mask, bool& carry) noexcept {
		const auto res = (val ^ mask) + carry;
		carry = res < carry;
		return res;
	}

	static constexpr auto Compare(uint64_t* a, uint64_t* b, size_t size) {
		if (std::is_constant_evaluated()) {
			return std::lexicographical_compare_three_way(a, a + size, b, b + size);
		} else {
			return std::memcmp(a, b, size * sizeof(uint64_t)) <=> 0;
		}
	}

public:
	static const Integer zero;
	static const Integer one;
	static const Integer minus_one;

private:
	// size_ == 0 => val_
	// size_ > 0 => bits_[size_]
	// size_ < 0 => bits_[-size_]
	ssize_t size_{};
	union {
		uint64_t* bits_;
		int64_t val_{};
	};

	[[nodiscard]] constexpr size_t size() const noexcept {
		return abs(size_);
	}

	[[nodiscard]] constexpr ssize_t ssize() const noexcept {
		return static_cast<ssize_t>(size());
	}

	[[nodiscard]] constexpr uint64_t* begin() const noexcept {
		return bits_;
	}

	[[nodiscard]] constexpr uint64_t* end() const noexcept {
		return bits_ + size();
	}

	[[nodiscard]] constexpr auto rbegin() const noexcept {
		return std::make_reverse_iterator(end());
	}

	[[nodiscard]] constexpr auto rend() const noexcept {
		return std::make_reverse_iterator(begin());
	}

	[[nodiscard]] constexpr uint64_t& front() const noexcept {
		return bits_[0];
	}

	[[nodiscard]] constexpr uint64_t& back() const noexcept {
		return bits_[size() - 1];
	}

	constexpr explicit Integer(const std::string& val, int radix = 10) {
		if (radix < 2 || radix > 36) {
			throw std::out_of_range("radix must be between 2 and 36");
		}

//		auto&& [sign_str, digits] = Check_number_string(val, radix);
//		if (digits == "0") {
//			size_ = 0;
//			val_ = 0;
//			return;
//		}
//		//TODO
//		const bool is_signed = sign_str == "-";
//		if (digits.size() <= digits_per_uint64[radix]) {
//			From(stoull(digits, nullptr, radix), is_signed);
//			return;
//		}
//
//		switch (radix) {
//		case 2:
//			From_radix2_4_16_string<2>(is_signed, digits);
//			return;
//		case 4:
//			From_radix2_4_16_string<4>(is_signed, digits);
//			return;
//		case 16:
//			From_radix2_4_16_string<16>(is_signed, digits);
//			return;
//		case 8:
//			From_radix8_32_string<8>(is_signed, digits);
//			return;
//		case 32:
//			From_radix8_32_string<32>(is_signed, digits);
//			return;
//		default:
//			From_string(is_signed, digits, radix);
//		}
	}

	constexpr Integer& operator=(const std::string& val) {
//		if (size_ != 0) {
//			delete[] bits_;
//		}
//
//		auto&& [sign_str, digits] = Check_number_string(val);
//		if (digits == "0") {
//			size_ = 0;
//			val_ = 0;
//			return *this;
//		}
//
//		const auto is_signed = sign_str == "-";
//		if (digits.size() <= digits_per_uint64[10]) {
//			From(stoull(digits), is_signed);
//			return *this;
//		}
//
//		From_string(is_signed, digits);
		return *this;
	}

	[[nodiscard]] std::string to_string(int radix = 10) const {
		if (radix < 2 || radix > 36) {
			throw std::out_of_range("radix must be between 2 and 36");
		}

		if (is_zero()) {
			return "0";
		}
//		//TODO
//		static constexpr char digits[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
//
//		const size_t n_bits = bit_count();
//		size_t n_digits = is_negative();
//		switch (radix) {
//		case 2:
//			n_digits += n_bits;
//			break;
//		case 4:
//			n_digits += (n_bits + 1) / 2;
//			break;
//		case 8:
//			n_digits += (n_bits + 2) / 3;
//			break;
//		case 16:
//			n_digits += (n_bits + 3) / 4;
//			break;
//		case 32:
//			n_digits += (n_bits + 4) / 5;
//			break;
//		default:
//			n_digits += static_cast<size_t>(ceil(static_cast<double>(n_bits) * digits_per_bit_factor[radix]));
//		}
//
//		std::string result;
//		result.reserve(n_digits);
//
//		if (size_ == 0) {
//			auto val = abs(val_);
//			do {
//				result += digits[val % radix];
//				val /= radix;
//			} while (val);
//		} else {
//			auto val = abs();
//			// TODO
//			do {
//				result += digits[static_cast<size_t>(val % radix)];
//				val /= radix;
//			} while (val);
//		}
//		if (is_negative()) {
//			result += '-';
//		}
//		std::reverse(result.begin(), result.end());
//		return result;
	}




public: // |
	friend constexpr Integer operator|(const Integer& lhs, const Integer& rhs) {
		if (lhs.is_zero()) {
			return rhs;
		}
		if (rhs.is_zero()) {
			return lhs;
		}

		if (lhs.size() <= 1 && rhs.size() <= 1) {
			return static_cast<__int128_t>(lhs) | static_cast<__int128_t>(rhs);
		}

		auto&& [longer, shorter] = Sort(lhs.ArrayStyle(), rhs.ArrayStyle());

		if (lhs.size_ > 1 && rhs.size_ > 1) {
			auto ret = Build(longer.size_);
			std::transform(
					longer.begin(), longer.begin() + shorter.size_, shorter.begin(), ret.begin(),
					std::bit_or<uint64_t>{}
			);
			std::copy(longer.begin() + shorter.size_, longer.begin() + longer.size_, ret.begin() + shorter.size_);
			return ret;
		}

		auto lc = longer.size_ < 0;
		auto sc = shorter.size_ < 0;
		bool rc = lc | sc;

		const uint64_t lm = -lc;
		const uint64_t sm = -sc;
		const uint64_t rm = -rc;

		const auto size = sc ? shorter.ssize() : longer.ssize();
		auto ret = Build(rc ? -size : size);
		std::transform(
				longer.begin(), longer.begin() + shorter.size(), shorter.begin(), ret.begin(),
				[&](auto l, auto s) { return Complement(Complement(l, lm, lc) | Complement(s, sm, sc), rm, rc); }
		);
		std::transform(
				longer.begin() + shorter.size(), longer.begin() + size, ret.begin() + shorter.size(),
				[&](auto l) { return Complement(Complement(l, lm, lc) | sm, rm, rc); }
		);

		if (rc) {
			ret.Extend1(1);
		} else {
			ret.Shrink();
		}
		return ret;
	}

	constexpr Integer& operator|=(const Integer& other) {
		if (this != &other || *this != other) {
			*this = *this | other;
		}
		return *this;
	}

public: // ^
	friend constexpr Integer operator^(const Integer& lhs, const Integer& rhs) {
		if (lhs.is_zero()) {
			return rhs;
		}
		if (rhs.is_zero()) {
			return lhs;
		}

		if (lhs.size() <= 1 && rhs.size() <= 1) {
			return static_cast<__int128_t>(lhs) ^ static_cast<__int128_t>(rhs);
		}

		auto&& [longer, shorter] = Sort(lhs.ArrayStyle(), rhs.ArrayStyle());

		if (lhs.size_ > 1 && rhs.size_ > 1) {
			if (lhs.size_ == rhs.size_) {
				auto size = lhs.size_;
				while (size > 0 && (lhs.bits_[size - 1] ^ rhs.bits_[size - 1]) == 0) {
					--size;
				}
				if (size == 0) {
					return zero;
				}
				if (size == 1) {
					return lhs.front() ^ rhs.front();
				}
				auto ret = Build(size);
				std::transform(lhs.begin(), lhs.end(), rhs.begin(), ret.begin(), std::bit_xor<uint64_t>{});
				return ret;
			}

			auto ret = Build(longer.size_);
			std::transform(
					longer.begin(), longer.begin() + shorter.size_, shorter.begin(), ret.begin(),
					std::bit_xor<uint64_t>{}
			);
			std::copy(longer.begin() + shorter.size_, longer.begin() + longer.size_, ret.begin() + shorter.size_);
			return ret;
		}

		auto lc = longer.size_ < 0;
		auto sc = shorter.size_ < 0;
		bool rc = lc ^ sc;

		const uint64_t lm = -lc;
		const uint64_t sm = -sc;
		const uint64_t rm = -rc;

		auto ret = Build(rc ? -longer.ssize() : longer.ssize());
		std::transform(
				longer.begin(), longer.begin() + shorter.size(), shorter.begin(), ret.begin(),
				[&](auto l, auto s) { return Complement(Complement(l, lm, lc) ^ Complement(s, sm, sc), rm, rc); }
		);
		std::transform(
				longer.begin() + shorter.size(), longer.begin() + longer.ssize(), ret.begin() + shorter.size(),
				[&](auto l) { return Complement(Complement(l, lm, lc), rm, rc); }
		);

		if (rc) {
			ret.Extend1(1);
		} else {
			ret.Shrink();
		}
		return ret;
	}

	constexpr Integer& operator^=(const Integer& other) {
		if (this == &other || *this == other) {
			*this = zero;
		} else {
			*this = *this ^ other;
		}
		return *this;
	}

public: // <<
	friend constexpr Integer operator<<(const Integer& lhs, size_t shift) {
		if (shift == 0 || lhs.is_zero()) {
			return lhs;
		}

		const auto q = shift / 64;
		const auto r = shift % 64;

		const auto size = static_cast<ssize_t>((lhs.AbsBitCount() + shift + 63) / 64);

		if (lhs.size_ == 0) {
			if (size == 1) {
				return { static_cast<uint64_t>(lhs.val_ << shift), lhs.val_ < 0 };
			}

			auto ret = Build(lhs.size_ < 0 ? -size : size);
			std::fill_n(ret.begin(), 0, size - 2);
			const auto val = abs(lhs.val_);
			ret.back() = val >> (64 - r);
			*(ret.end() - 2) = val << r;
			return ret;
		}

		auto ret = Build(lhs.size_ < 0 ? -size : size);
		std::fill_n(ret.begin(), 0, q);
		if (r == 0) {
			std::copy(lhs.begin() + q, lhs.end(), ret.begin() + q);
		} else {
			uint64_t carry = 0;
			std::transform(
					lhs.begin() + q, lhs.end(), ret.begin() + q, [&](auto e) {
						carry = e >> (64 - r);
						return e << r | carry;
					}
			);
			if (carry) {
				ret.back() = carry;
			}
		}
		return ret;
	}

	constexpr Integer& operator<<=(size_t shift) {
		*this = *this << shift;
		return *this;
	}

public: // >>
	friend constexpr Integer operator>>(const Integer& lhs, size_t shift) {
		if (shift == 0 || lhs.is_zero()) {
			return lhs;
		}

		if (lhs.size_ == 0) {
			return lhs.val_ >> shift;
		}

		const auto bit_count = lhs.AbsBitCount();
		if (shift >= bit_count) {
			return lhs.size_ > 0 ? zero : minus_one;
		}

		const auto q = shift / 64;
		const auto r = shift % 64;

		const auto size = (static_cast<ssize_t>(bit_count - shift + 63) / 64);

		if (size == 1) {
			if (lhs.size() == 1) {
				assert(static_cast<int64_t>(lhs.front()) < 0);
				const auto val = static_cast<int64_t>(lhs.front() >> shift);
				return lhs.size_ < 0 ? -val : val;
			}

			const auto val = lhs.back() << (64 - r) | *(lhs.end() - 2) >> r;
			return { val, lhs.size_ < 0 };
		}

		auto ret = Build(lhs.size_ < 0 ? -size : size);
		if (r == 0) {
			std::copy(lhs.begin() + q, lhs.end(), ret.begin());
		} else {
			uint64_t carry = lhs.front() >> r;
			std::transform(
					lhs.begin() + 1, lhs.end(), ret.begin(), [&](auto e) {
						carry = e >> r;
						return e << (64 - r) | carry;
					}
			);
			ret.back() = carry;
		}
		return ret;
	}

	constexpr Integer& operator>>=(size_t shift) {
		*this = *this >> shift;
		return *this;
	}

public: // +
	friend constexpr Integer operator+(const Integer& lhs, const Integer& rhs) {
		if (lhs.is_zero()) {
			return rhs;
		}
		if (rhs.is_zero()) {
			return lhs;
		}
		return lhs.sign() == rhs.sign() ? Add(lhs, rhs) : Sub(lhs, rhs.RefNeg());
	}

	constexpr Integer& operator+=(const Integer& other) {
		*this = *this + other;
		return *this;
	}

public: // -
	friend constexpr Integer operator-(const Integer& lhs, const Integer& rhs) {
		if (lhs.is_zero()) {
			return -rhs;
		}
		if (rhs.is_zero()) {
			return lhs;
		}
		return lhs.sign() == rhs.sign() ? Sub(lhs, rhs) : Add(lhs, rhs.RefNeg());
	}

	constexpr Integer& operator-=(const Integer& other) {
		*this = *this - other;
		return *this;
	}

public: // *
	friend constexpr Integer operator*(const Integer& lhs, const Integer& rhs) {
		if (lhs.is_zero() || rhs.is_zero()) {
			return zero;
		}

		bool negative = (lhs.sign() ^ rhs.sign()) < 0;

		if (lhs.size() <= 1 && rhs.size() <= 1) {
			return { lhs.U128Abs() * rhs.U128Abs(), negative };
		}

		if (lhs.size() <= 1) {
			return Mul1(rhs, lhs.U64Abs(), negative);
		}
		if (rhs.size() <= 1) {
			return Mul1(lhs, rhs.U64Abs(), negative);
		}

		return lhs.size() > rhs.size() ? Mul2(lhs, rhs, negative) : Mul2(rhs, lhs, negative);
	}

	constexpr Integer& operator*=(const Integer& other) {
		//TODO
		return *this;
	}

public: // /
	friend constexpr Integer operator/(const Integer& dividend, const Integer& divisor) {
		if (divisor.is_zero()) {
			throw std::invalid_argument("division by zero");
		}

		const auto cmp = dividend.RefAbs() <=> divisor.RefAbs();

		if (cmp < nullptr) {
			return zero;
		}

		bool negative = (dividend.sign() ^ divisor.sign()) < 0;

		if (cmp == nullptr) {
			return negative ? minus_one : one;
		}

		if (dividend.size() <= 1) { //divisor.size() <= 1
			return { dividend.U64Abs() / divisor.U64Abs(), negative };
		}
		if (dividend.size() == 2) { //divisor.size() <= 2
			return { dividend.U128Abs() / divisor.U128Abs(), negative };
		}
		if (divisor.size() <= 1) { //dividend.size() > 2
			return Div1(dividend, divisor.U64Abs(), negative);
		}
		//dividend.size() > 2 divisor.size() >= 2
		return DivMod(dividend, divisor, negative).first;
	}

	constexpr Integer& operator/=(const Integer& other) {
		*this = *this / other;
		return *this;
	}

public: // %
	friend constexpr Integer operator%(const Integer& dividend, const Integer& divisor) {
		if (divisor.is_zero()) {
			throw std::invalid_argument("division by zero");
		}

		const auto cmp = dividend.RefAbs() <=> divisor.RefAbs();

		if (cmp < nullptr) {
			return dividend;
		}
		if (cmp == nullptr) {
			return zero;
		}

		bool negative = dividend.sign() < 0;

		if (dividend.size() <= 1) { //divisor.size() <= 1
			return { dividend.U64Abs() % divisor.U64Abs(), negative };
		}
		if (dividend.size() == 2) { //divisor.size() <= 2
			return { dividend.U128Abs() % divisor.U128Abs(), negative };
		}
		if (divisor.size() <= 1) { //dividend.size() > 2
			return Mod1(dividend, divisor.U64Abs(), negative);
		}
		return DivMod(dividend, divisor, negative).second;
	}

	constexpr Integer& operator%=(const Integer& other) {
		*this = *this % other;
		return *this;
	}

public:
	friend std::istream& operator>>(std::istream& is, Integer& val) {
		std::string str;
		std::getline(is, str);
		val = str;
		return is;
	}

	friend std::ostream& operator<<(std::ostream& os, const Integer& val) {
		return os << val.to_string();
	}

private:
	constexpr void From(uint64_t val) {
		if (val <= std::numeric_limits<int64_t>::max()) {
			size_ = 0;
			val_ = static_cast<int64_t>(val);
		} else {
			size_ = 1;
			bits_ = new uint64_t{ val };
		}
	}

	constexpr void From(__uint128_t val) {
		if (val <= std::numeric_limits<int64_t>::max()) {
			size_ = 0;
			val_ = static_cast<int64_t>(val);
		} else if (val <= std::numeric_limits<uint64_t>::max()) {
			size_ = 1;
			bits_ = new uint64_t{ static_cast<uint64_t>(val) };
		} else {
			size_ = 2;
			bits_ = new uint64_t[2]{ static_cast<uint64_t>(val), static_cast<uint64_t>(val >> 64) };
		}
	}

	constexpr void Extend1(uint64_t val) {
		size_ = size_ < 0 ? size_ - 1 : size_ + 1;
		bits_ = static_cast<uint64_t*>(std::realloc(bits_, size()));
		back() = val;
	}

	constexpr void Shrink() {
		const auto p = std::find_if_not(rbegin(), rend(), [](auto e) { return e == 0; });
		const auto size = rend() - p;
		if (size == 0) {
			*this = zero;
		} else if (size == 1) {
			*this = Integer(front(), size_ < 0);
		} else if (size < this->size()) {
			auto ret = Build(size_ < 0 ? -size : size);
			std::copy_n(begin(), size, ret.begin());
			*this = std::move(ret);
		}
	}

	[[nodiscard]] constexpr uint64_t U64Abs() const {
		return size_ == 0 ? abs(val_) : front();
	}

	[[nodiscard]] constexpr __uint128_t U128Abs() const {
		if (size() <= 1) {
			return static_cast<__uint128_t>(U64Abs());
		}
		return static_cast<__uint128_t>(*(begin() + 1)) << 64 | front();
	}

	[[nodiscard]] constexpr Integer ArrayStyle() const {
		if (size_ != 0) {
			return *this;
		}
		return { new uint64_t{ static_cast<uint64_t>(abs(val_)) }, val_ < 0 ? -1 : 1 };
	}

	[[nodiscard]] constexpr Integer RefNeg() const {
		if (size_ == 0) {
			return { static_cast<uint64_t>(-val_), val_ > 0 };
		}
		return { bits_, -size_ };
	}

	[[nodiscard]] constexpr Integer RefAbs() const {
		if (size_ == 0) {
			return abs(val_);
		}
		return { bits_, ssize() };
	}

	[[nodiscard]] constexpr size_t AbsBitCount() const {
		if (size_ == 0) {
			return bit_count(abs(val_));
		}
		return (size() - 1) * 64 + bit_count(back());
	}

private:
	static constexpr std::pair<const Integer&, const Integer&> Sort(const Integer& a, const Integer& b) noexcept {
		return a.size() > b.size() ? std::make_pair(a, b) : std::make_pair(b, a);
	}

	// a.sign == b.sign
	static constexpr Integer Add(const Integer& a, const Integer& b) {
		if (a.size() <= 1 && b.size() <= 1) {
			return static_cast<__int128_t>(a) + static_cast<__int128_t>(b);
		}

		bool negative = a.sign() < 0;

		if (a.size() <= 1) {
			return Add1(b, a.U64Abs(), negative);
		}
		if (b.size() <= 1) {
			return Add1(a, b.U64Abs(), negative);
		}

		return a.size() > b.size() ? Add2(a, b, negative) : Add2(b, a, negative);
	}

	static constexpr Integer Add1(const Integer& a, uint64_t b, bool negative) {
		auto ret = Build(negative ? -a.ssize() : a.ssize());
		bool carry = false;
		ret.front() = Add(a.front(), b, carry);
		std::transform(
				a.begin() + 1, a.end(), ret.begin() + 1, [&carry](auto e) {
					return Add(e, carry);
				}
		);
		if (carry) {
			ret.Extend1(1);
		}
		return ret;
	}

	static constexpr Integer Add2(const Integer& a, const Integer& b, bool negative) {
		auto ret = Build(negative ? -a.ssize() : a.ssize());
		bool carry = false;
		std::transform(
				a.begin(), a.begin() + b.size(), b.begin(), ret.begin(), [&carry](auto l, auto r) {
					return Add(l, r, carry);
				}
		);
		std::transform(
				a.begin() + b.size(), a.end(), ret.begin() + b.size(), [&carry](auto e) {
					return Add(e, carry);
				}
		);
		if (carry) {
			ret.Extend1(1);
		}
		return ret;
	}

	// a.sign == b.sign
	static constexpr Integer Sub(const Integer& a, const Integer& b) {
		const auto cmp = a.RefAbs() <=> b.RefAbs();

		if (cmp == nullptr) {
			return zero;
		}

		if (a.size() <= 1 && b.size() <= 1) {
			return static_cast<__int128_t>(a) - static_cast<__int128_t>(b);
		}

		if (a.size() <= 1) {
			return Sub1(b, a.U64Abs(), a.sign() > 0);
		}
		if (b.size() <= 1) {
			return Sub1(a, b.U64Abs(), a.sign() < 0);
		}

		return cmp > nullptr ? Sub2(a, b, a.size_ < 0) : Sub2(b, a, a.size_ > 0);
	}

	static constexpr Integer Sub1(const Integer& a, uint64_t b, bool negative) {
		auto ret = Build(negative ? -a.ssize() : a.ssize());
		bool borrow = false;
		ret.front() = Sub(a.front(), b, borrow);
		std::transform(
				a.begin() + 1, a.end(), ret.begin() + 1, [&borrow](auto e) {
					return Sub(e, borrow);
				}
		);
		ret.Shrink();
		return ret;
	}

	static constexpr Integer Sub2(const Integer& a, const Integer& b, bool negative) {
		auto ret = Build(negative ? -a.ssize() : a.ssize());
		bool borrow = false;
		std::transform(
				a.begin(), a.begin() + b.size(), b.begin(), ret.begin(), [&borrow](auto l, auto r) {
					return Sub(l, r, borrow);
				}
		);
		std::transform(
				a.begin() + b.size(), a.end(), ret.begin() + b.size(), [&borrow](auto e) {
					return Sub(e, borrow);
				}
		);
		ret.Shrink();
		return ret;
	}

	static constexpr Integer Mul1(const Integer& a, uint64_t b, bool negative) {
		auto ret = Build(negative ? -a.ssize() : a.ssize());
		uint64_t carry = 0;
		std::transform(
				a.begin(), a.end(), ret.begin(), [&](auto e) {
					return MulAdd(e, b, carry);
				}
		);
		if (carry) {
			ret.Extend1(carry);
		}
		return ret;
	}

	static constexpr Integer Mul2(const Integer& a, const Integer& b, bool negative) {
		auto ret = Build(negative ? -(a.ssize() + b.ssize()) : a.ssize() + b.ssize());

		uint64_t carry = 0;
		*std::transform(
				a.begin(), a.end(), ret.begin(), [&](auto e) { return MulAdd(e, ret.front(), carry); }
		) = carry;

		for (auto p = b.begin() + 1, r = ret.begin() + 1; p != b.end(); ++p, ++r) {
			const auto y = *p;
			carry = 0;
			*std::transform(
					a.begin(), a.end(), r, r, [&](auto x, auto last_sum) { return MulAdd(x, y, last_sum, carry); }
			) = carry;
		}

		ret.Shrink();
		return ret;
	}

	static constexpr Integer Div1(const Integer& dividend, uint64_t divisor, bool negative) {
		const auto size = dividend.back() >= divisor ? dividend.ssize() : dividend.ssize() - 1;
		auto ret = Build(negative ? -size : size);
		uint64_t high = 0;
		std::transform(
				dividend.rbegin(), dividend.rend(), ret.rbegin(), [&](auto low) {
					return Div(high, low, divisor);
				}
		);
		return ret;
	}

	static constexpr Integer Mod1(const Integer& dividend, uint64_t divisor, bool negative) {
		uint64_t high = 0;
		std::for_each(
				dividend.rbegin(), dividend.rend(), [&](auto low) {
					high = (static_cast<__uint128_t>(high) << 64 | low) % divisor;
				}
		);
		return high;
	}

	static constexpr std::pair<Integer, Integer> DivMod(const Integer& dividend, const Integer& divisor, bool negative) {
		return {};
	}


private:
	static constexpr double bits_per_digit_factor[] = {
			0, 0,
			1.00000000000000, 1.58496250072116, 2.00000000000000, 2.32192809488736, 2.58496250072116,
			2.80735492205760, 3.00000000000000, 3.16992500144231, 3.32192809488736, 3.45943161863730,
			3.58496250072116, 3.70043971814109, 3.80735492205760, 3.90689059560852, 4.00000000000000,
			4.08746284125034, 4.16992500144231, 4.24792751344359, 4.32192809488736, 4.39231742277876,
			4.45943161863730, 4.52356195605701, 4.58496250072116, 4.64385618977472, 4.70043971814109,
			4.75488750216347, 4.80735492205760, 4.85798099512757, 4.90689059560852, 4.95419631038687,
			5.00000000000000, 5.04439411935845, 5.08746284125034, 5.12928301694497, 5.16992500144231
	};

	static constexpr double digits_per_bit_factor[] = {
			0, 0,
			1.000000000000000, 0.630929753571457, 0.500000000000000, 0.430676558073393, 0.386852807234542,
			0.356207187108022, 0.333333333333333, 0.315464876785729, 0.301029995663981, 0.289064826317888,
			0.278942945651130, 0.270238154427320, 0.262649535037194, 0.255958024809815, 0.250000000000000,
			0.244650542118226, 0.239812466568131, 0.235408913366638, 0.231378213159759, 0.227670248696953,
			0.224243824217575, 0.221064729457504, 0.218104291985532, 0.215338279036697, 0.212746053553363,
			0.210309917857152, 0.208014597676509, 0.205846832460434, 0.203795047090506, 0.201849086582100,
			0.200000000000000, 0.198239863170561, 0.196561632232823, 0.194959021893786, 0.193426403617271
	};

	static constexpr size_t digits_per_uint64[] = {
			0, 0,
			64, 40, 32, 27, 24,
			22, 21, 20, 19, 18,
			17, 17, 16, 16, 16,
			15, 15, 15, 14, 14,
			14, 14, 13, 13, 13,
			13, 13, 13, 13, 12,
			12, 12, 12, 12, 12
	};

	static constexpr size_t super_radix[] = {
			0, 0,
			0xffffffffffffffff, 0xa8b8b452291fe820, 0xffffffffffffffff, 0x6765c793fa10079c, 0x41c21cb8e0ffffff,
			0x3642798750226110, 0x7fffffffffffffff, 0xa8b8b452291fe820, 0x8ac7230489e7ffff, 0x4d28cb56c33fa538,
			0x1eca170bffffffff, 0x780c7372621bd74c, 0x1e39a5057d80ffff, 0x5b27ac993df97700, 0xffffffffffffffff,
			0x27b95e997e21d9f0, 0x5da0e1e53c5c7fff, 0xd2ae3299c1c4aeda, 0x16bcc41e8fffffff, 0x2d04b7fdd9c0ef48,
			0x5658597bcaa23fff, 0xa0e2073737609370, 0x0c29e97fffffffff, 0x14adf4b7320334b8, 0x226ed36478bf9fff,
			0x383d9170b85ff80a, 0x5a3c23e39bffffff, 0x8e65137388122bcc, 0xdd41bb36d259dfff, 0x0aee5720ee830680,
			0x0fffffffffffffff, 0x172588ad4f5f0980, 0x211e44f7d02c0fff, 0x2ee56725f06e5c70, 0x41c21cb8e0ffffff
	};


//	template <size_t radix>
//	constexpr void From_radix2_4_16_string(bool is_signed, const std::string& digits) {
//		const ssize_t len = ((digits.size() * tailing_zeros(radix) + 63) & ~63) / 64;
//		size_ = is_signed ? -len : len;
//		bits_ = New(len);
//		constexpr auto step = digits_per_uint64[radix];
//		auto p = bits_;
//		for (auto it = digits.end(); it > digits.begin(); it -= step) {
//			*(p++) = stoull(std::string(max(digits.begin(), it - step), it), nullptr, radix);
//		}
//	}
//
//	template <size_t radix>
//	constexpr void From_radix8_32_string(bool is_signed, const std::string& digits) {
//		const ssize_t prob_len = ((digits.size() * tailing_zeros(radix) + 63) & ~63) / 64 + 1;
//		const auto bits = New(prob_len);
//		auto last = bits + prob_len;
//		constexpr auto step = digits_per_uint64[radix];
//
//		auto p = bits;
//		auto it = digits.end();
//		*(p++) = stoull(std::string(it - step, it), nullptr, radix);
//
//		constexpr auto diff = 64 - step * tailing_zeros(radix);
//		auto shift = diff;
//		for (it -= step; it > digits.begin(); it -= step) {
//			auto num = stoull(std::string(max(digits.begin(), it - step), it), nullptr, radix);
//			if (shift != 8) {
//				*(p - 1) |= num << (64 - shift);
//				*(p++) = num >> shift;
//				shift += diff;
//			} else {
//				p = reinterpret_cast<uint64_t*>(reinterpret_cast<uint8_t*>(p) + 7);
//				*(p++) = num;
//				shift = diff;
//			}
//		}
//		std::fill(reinterpret_cast<uint8_t*>(p), reinterpret_cast<uint8_t*>(last), 0);
//
//		while (*(last - 1) == 0) {
//			--last;
//		}
//		ssize_t len = last - bits;
//		size_ = is_signed ? -len : len;
//		bits_ = New(len);
//		std::copy_n(bits, len, bits_);
//		delete[] bits;
//	}
//
//	void From_string(bool is_signed, const std::string& digits, int radix = 10) {
//		const auto n_bits = static_cast<size_t>(ceil(
//				static_cast<double>(digits.size()) * bits_per_digit_factor[radix]
//		));
//		const ssize_t prob_len = ((n_bits + 63) & ~63) / 64;
//		const auto first = New(prob_len);
//		auto last = first + prob_len;
//		const auto step = digits_per_uint64[radix];
//
//		auto p = reinterpret_cast<unsigned long long*>(first);
//		auto it = digits.begin();
//		*(p++) = stoull(std::string(it, it + step), nullptr, radix);
//
//		for (it += step; it < digits.end(); it += step) {
//			auto carry = stoull(std::string(it, it + step), nullptr, radix);
//			const auto q = std::exchange(p, reinterpret_cast<unsigned long long*>(first));
//			while (p != q) {
//				MulAdd(p++, super_radix[radix], carry);
//			}
//			if (carry != 0) {
//				*(p++) = carry;
//			}
//		}
//		std::fill(reinterpret_cast<uint64_t*>(p), last, 0);
//
//		while (*(last - 1) == 0) {
//			--last;
//		}
//		ssize_t len = last - first;
//		if (len == 1) {
//			From(is_signed, *first);
//			delete[] first;
//		} else if (len == prob_len) {
//			size_ = is_signed ? -prob_len : prob_len;
//			bits_ = first;
//		} else {
//			size_ = is_signed ? -len : len;
//			bits_ = New(len);
//			std::copy_n(first, len, bits_);
//			delete[] first;
//		}
//	}


//	static std::pair<std::string, std::string> Check_number_string(const std::string& val, int radix = 10) {
//		if (val.empty()) {
//			throw std::invalid_argument("zero length string");
//		}
//
//		std::stringstream ss;
//		ss << "[ \t]*([+-])?0*([0-";
//		if (radix <= 10) {
//			ss << radix - 1;
//		} else {
//			char d = static_cast<char>(radix - 11);
//			ss << "9a-" << 'a' + d << "A-" << 'A' + d;
//		}
//		ss << "]+)[ \t]*";
//
//		std::smatch matches;
//		if (!std::regex_match(val, matches, std::regex(ss.str()))) {
//			throw std::invalid_argument("illegal digit");
//		}
//
//		return { matches[1].str(), matches[2].str() };
//	}
};


const Integer Integer::zero = Integer();
const Integer Integer::one = Integer(1);
const Integer Integer::minus_one = Integer(-1);

