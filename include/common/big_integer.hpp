
//
// (C) Daniel Strano and the Qimcifa contributors, 2022, 2023. All rights reserved.
// (C) Quantum Circuits, Inc, 2024. All rights reserved.
//
// This header has been adapted for OpenCL and C, from big_integer.c by Andre Azevedo.
// This header has been restored to sanity by Stefan Teleman for Quantum Circuits, Inc.
// It now uses GNU MP for Big Integers.
//
// Original file:
//
// big_integer.c
//     Description: "Arbitrary"-precision integer
//     Author: Andre Azevedo <http://github.com/andreazevedo>
//
// The MIT License (MIT)
//
// Copyright (c) 2014 Andre Azevedo
// Copyright (c) 2024 Quantum Circuits, Inc.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <climits>
#include <cassert>
#include <type_traits>
#include <stdexcept>

#include <gmp.h>

namespace Qrack {

#define BIG_INTEGER_WORD_BITS 64U
#define BIG_INTEGER_WORD_POWER 6U
#define BIG_INTEGER_WORD uint64_t
#define BIG_INTEGER_HALF_WORD uint32_t
#define BIG_INTEGER_HALF_WORD_POW 0x100000000ULL
#define BIG_INTEGER_HALF_WORD_MASK 0xFFFFFFFFULL
#define BIG_INTEGER_HALF_WORD_MASK_NOT 0xFFFFFFFF00000000ULL

// This can be any power of 2 greater than (or equal to) 64:
#define BIG_INTEGER_BITS (1U << QBCAPPOW)
#define BIG_INTEGER_WORD_SIZE (int64_t)(BIG_INTEGER_BITS / BIG_INTEGER_WORD_BITS)

// The rest of the constants need to be consistent with the one above:
constexpr size_t BIG_INTEGER_HALF_WORD_BITS = BIG_INTEGER_WORD_BITS >> 1U;
constexpr int32_t BIG_INTEGER_HALF_WORD_SIZE = BIG_INTEGER_WORD_SIZE << 1U;
constexpr int32_t BIG_INTEGER_MAX_WORD_INDEX = BIG_INTEGER_WORD_SIZE - 1U;
constexpr int32_t MPZ_INTEGER_BITS = BIG_INTEGER_WORD_SIZE << UINTPOW;

class BigInteger {
public:
  mutable mpz_t mpz;

public:
  BigInteger() : mpz() {
    mpz_init2(mpz, MPZ_INTEGER_BITS);
    mpz_set_ui(mpz, 0UL);
  }

  BigInteger(const BigInteger& rhs) : mpz() {
    mpz_init2(mpz, MPZ_INTEGER_BITS);
    mpz_set(mpz, rhs.mpz);
  }

  BigInteger(BigInteger&& rhs) : mpz() {
    mpz_init2(mpz, MPZ_INTEGER_BITS);
    mpz_set(mpz, rhs.mpz);
    assert(mpz_cmp(mpz, rhs.mpz) == 0 && "invalid mpz!");
  }

  BigInteger(BIG_INTEGER_WORD rhs) : mpz() {
    mpz_init2(mpz, MPZ_INTEGER_BITS);
    mpz_set_ui(mpz, static_cast<uint64_t>(rhs));
  }

  BigInteger(const mpz_t& rhs) : mpz() {
    mpz_init2(mpz, MPZ_INTEGER_BITS);
    mpz_set(mpz, rhs);
    assert(mpz_cmp(mpz, rhs) == 0 && "invalid mpz!");
  }

  BigInteger(mpz_t&& rhs) : mpz() {
    mpz_init2(mpz, MPZ_INTEGER_BITS);
    mpz_set(mpz, rhs);
    assert(mpz_cmp(mpz, rhs) == 0 && "invalid mpz!");
  }

  ~BigInteger() {
    (void) mpz_clear(mpz);
  }

  BigInteger& operator=(const BigInteger& rhs) {
    if (this != &rhs) {
      mpz_set(mpz, rhs.mpz);
      assert(mpz_cmp(mpz, rhs.mpz) == 0 && "invalid mpz!");
    }

    return *this;
  }

  BigInteger& operator=(BigInteger&& rhs) {
    if (this != &rhs) {
      mpz_set(mpz, rhs.mpz);
      assert(mpz_cmp(mpz, rhs.mpz) == 0 && "invalid mpz!");
    }

    return *this;
  }

  BigInteger& operator=(BIG_INTEGER_WORD rhs) {
    mpz_set_ui(mpz, rhs);
    return *this;
  }

  BigInteger& operator=(const mpz_t& rhs) {
    mpz_set(mpz, rhs);
    assert(mpz_cmp(mpz, rhs) == 0 && "invalid mpz!");
    return *this;
  }

  explicit inline operator BIG_INTEGER_WORD() {
    return mpz_get_ui(mpz);
  }

  explicit inline operator BIG_INTEGER_WORD() const {
    return mpz_get_ui(mpz);
  }

  explicit inline operator bool() {
    return mpz_cmp_ui(mpz, 0UL) != 0;
  }

  explicit inline operator bool() const {
    return mpz_cmp_ui(mpz, 0UL) != 0;
  }

  explicit inline operator uint32_t() {
    return static_cast<uint32_t>(mpz_get_ui(mpz));
  }

  explicit inline operator uint32_t() const {
    return static_cast<uint32_t>(mpz_get_ui(mpz));
  }

  explicit inline operator uint16_t() {
    return static_cast<uint16_t>(mpz_get_ui(mpz));
  }

  explicit inline operator uint16_t() const {
    return static_cast<uint16_t>(mpz_get_ui(mpz));
  }

  explicit inline operator uint8_t() {
    return static_cast<uint8_t>(mpz_get_ui(mpz));
  }

  explicit inline operator uint8_t() const {
    return static_cast<uint8_t>(mpz_get_ui(mpz));
  }

  explicit inline operator double() {
    return mpz_get_d(mpz);
  }

  explicit inline operator double() const {
    return mpz_get_d(mpz);
  }

  // Standard Arithmetic Operators.
  inline BigInteger operator+(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_add(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator+(BIG_INTEGER_WORD rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_add_ui(result, mpz, rhs);
    return result;
  }

  inline BigInteger operator+(BIG_INTEGER_HALF_WORD rhs) const {
    return BigInteger::operator+(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger operator-(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_sub(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator-(BIG_INTEGER_WORD rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_sub_ui(result, mpz, rhs);
    return result;
  }

  inline BigInteger operator-(BIG_INTEGER_HALF_WORD rhs) const {
    return BigInteger::operator-(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger operator*(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_mul(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator*(BIG_INTEGER_WORD rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_mul_ui(result, mpz, rhs);
    return result;
  }

  inline BigInteger operator*(BIG_INTEGER_HALF_WORD rhs) const {
    return BigInteger::operator*(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger operator/(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);

    if (mpz_sgn(mpz) == -1 && mpz_sgn(rhs.mpz) == -1) {
      mpz_cdiv_q(result, mpz, rhs.mpz);
    } else if (mpz_sgn(mpz) == -1 || mpz_sgn(rhs.mpz) == -1) {
      mpz_fdiv_q(result, mpz, rhs.mpz);
    } else if (mpz_sgn(mpz) == 0 || mpz_sgn(rhs.mpz) == 0) {
      return result;
    } else {
      mpz_cdiv_q(result, mpz, rhs.mpz);
    }

    return result;
  }

  inline BigInteger operator/(BIG_INTEGER_WORD rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);

    if (mpz_sgn(mpz) == -1) {
      (void) mpz_fdiv_q_ui(result, mpz, rhs);
    } else if (mpz_sgn(mpz) == 0 || rhs == 0) {
      return result;
    } else {
      (void) mpz_cdiv_q_ui(result, mpz, rhs);
    }

    return result;
  }

  inline BigInteger operator/(BIG_INTEGER_HALF_WORD rhs) const {
    return BigInteger::operator/(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger operator%(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);

    if (mpz_sgn(mpz) == -1 && mpz_sgn(rhs.mpz) == -1) {
      (void) mpz_cdiv_r(result, mpz, rhs.mpz);
    } else if (mpz_sgn(mpz) == -1 || mpz_sgn(rhs.mpz) == -1) {
      (void) mpz_fdiv_r(result, mpz, rhs.mpz);
    } else if (mpz_sgn(mpz) == 0 || mpz_sgn(rhs.mpz) == 0) {
      return result;
    } else {
      (void) mpz_cdiv_r(result, mpz, rhs.mpz);
    }

    return result;
  }

  inline BigInteger operator%(BIG_INTEGER_WORD rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);

    if (mpz_sgn(mpz) == -1) {
      (void) mpz_fdiv_r_ui(result, mpz, rhs);
    } else if (mpz_sgn(mpz) == 0 || rhs == 0) {
      return result;
    } else {
      (void) mpz_cdiv_r_ui(result, mpz, rhs);
    }

    return result;
  }

  inline BigInteger operator%(BIG_INTEGER_HALF_WORD rhs) const {
    return BigInteger::operator%(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger operator~() const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_mul_si(result, mpz, -1L);
    return result;
  }

  inline bool operator!() {
    return mpz_sgn(mpz) == 0;
  }

  inline bool operator==(const BigInteger& rhs) const {
    return mpz_cmp(mpz, rhs.mpz) == 0;
  }

  inline bool operator==(BIG_INTEGER_WORD rhs) const {
    return mpz_cmp_ui(mpz, rhs) == 0;
  }

  inline bool operator==(BIG_INTEGER_HALF_WORD rhs) const {
    return mpz_cmp_ui(mpz, rhs) == 0;
  }

  inline bool operator!=(const BigInteger& rhs) const {
    return mpz_cmp(mpz, rhs.mpz) != 0;
  }

  inline bool operator!=(BIG_INTEGER_WORD rhs) const {
    return mpz_cmp_ui(mpz, rhs) != 0;
  }

  inline bool operator!=(BIG_INTEGER_HALF_WORD rhs) const {
    return mpz_cmp_ui(mpz, rhs) != 0;
  }

  inline bool operator<(const BigInteger& rhs) const {
    return mpz_cmp(mpz, rhs.mpz) < 0;
  }

  inline bool operator<(BIG_INTEGER_WORD rhs) const {
    return mpz_cmp_ui(mpz, rhs) < 0;
  }

  inline bool operator>(const BigInteger& rhs) const {
    return mpz_cmp(mpz, rhs.mpz) > 0;
  }

  inline bool operator>(BIG_INTEGER_WORD rhs) const {
    return mpz_cmp_ui(mpz, rhs) > 0;
  }

  template<typename __ST>
  inline BigInteger operator<<(__ST sq) const {
    static_assert(std::is_integral<__ST>::value,
                  "shift quantity must be an integral type!");
    static_assert(!std::is_signed<__ST>::value,
                  "shift quantity must be an unsigned integral type!");
    if (!sq)
      return BigInteger(*this);

    if (this->is_zero())
      return BigInteger(0UL);

    mpz_t r;
    mpz_init2(r, MPZ_INTEGER_BITS);

    for (__ST i = 0; i < sq + 1U; ++i) {
      mpz_mul_2exp(r, mpz, i);
    }

    return BigInteger(r);
  }

  template<typename __ST>
  inline BigInteger operator>>(__ST sq) const {
    static_assert(std::is_integral<__ST>::value,
                  "shift quantity must be an integral type!");
    static_assert(!std::is_signed<__ST>::value,
                  "shift quantity must be an unsigned integral type!");
    if (!sq)
      return BigInteger(*this);

    if (this->is_zero())
      return BigInteger(0UL);

    mpz_t r;
    mpz_init2(r, MPZ_INTEGER_BITS);

    if (this->is_negative()) {
      for (__ST i = 0; i < sq + 1U; ++i) {
        mpz_mul_2exp(r, mpz, i);
        mpz_fdiv_q_2exp(r, mpz, i);
      }
    } else {
      for (__ST i = 0; i < sq + 1U; ++i) {
        mpz_mul_2exp(r, mpz, i);
        mpz_tdiv_q_2exp(r, mpz, i);
      }
    }

    return BigInteger(r);
  }

  inline BigInteger operator&(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_and(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator^(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_xor(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator|(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_ior(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator&(BIG_INTEGER_WORD rhs) const {
    mpz_t op;
    mpz_t result;
    mpz_init2(op, MPZ_INTEGER_BITS);
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_set_ui(op, rhs);
    mpz_and(result, mpz, op);
    mpz_clear(op);
    return result;
  }

  inline BigInteger operator^(BIG_INTEGER_WORD rhs) const {
    mpz_t op;
    mpz_t result;
    mpz_init2(op, MPZ_INTEGER_BITS);
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_set_ui(op, rhs);
    mpz_xor(result, mpz, op);
    mpz_clear(op);
    return result;
  }

  inline BigInteger operator|(BIG_INTEGER_WORD rhs) const {
    mpz_t op;
    mpz_t result;
    mpz_init2(op, MPZ_INTEGER_BITS);
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_set_ui(op, rhs);
    mpz_ior(result, mpz, op);
    mpz_clear(op);
    return result;
  }

  inline BigInteger& operator++() {
    mpz_add_ui(mpz, mpz, 1UL);
    return *this;
  }

  inline BigInteger& operator--() {
    mpz_sub_ui(mpz, mpz, 1UL);
    return *this;
  }

  inline BigInteger operator++(int) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_set(result, mpz);
    mpz_add_ui(mpz, mpz, 1UL);
    return result;
  }

  inline BigInteger operator--(int) const {
    mpz_t result;
    mpz_init2(result, MPZ_INTEGER_BITS);
    mpz_set(result, mpz);
    mpz_sub_ui(mpz, mpz, 1UL);
    return result;
  }

  inline BigInteger& operator+=(const BigInteger& rhs) {
    mpz_add(mpz, rhs.mpz, mpz);
    return *this;
  }

  inline BigInteger& operator+=(BIG_INTEGER_WORD rhs) {
    mpz_add_ui(mpz, mpz, rhs);
    return *this;
  }

  inline BigInteger& operator+=(BIG_INTEGER_HALF_WORD rhs) {
    return BigInteger::operator+=(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger& operator-=(const BigInteger& rhs) {
    mpz_sub(mpz, mpz, rhs.mpz);
    return *this;
  }

  inline BigInteger& operator-=(BIG_INTEGER_WORD rhs) {
    mpz_sub_ui(mpz, mpz, rhs);
    return *this;
  }

  inline BigInteger& operator-=(BIG_INTEGER_HALF_WORD rhs) {
    return BigInteger::operator-=(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger& operator*=(const BigInteger& rhs) {
    mpz_mul(mpz, mpz, rhs.mpz);
    return *this;
  }

  inline BigInteger& operator*=(BIG_INTEGER_WORD rhs) {
    mpz_mul_ui(mpz, mpz, rhs);
    return *this;
  }

  inline BigInteger& operator*=(BIG_INTEGER_HALF_WORD rhs) {
    return BigInteger::operator*=(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger& operator/=(const BigInteger& rhs) {
    if (mpz_sgn(mpz) == -1 && mpz_sgn(rhs.mpz) == -1) {
      mpz_cdiv_q(mpz, mpz, rhs.mpz);
    } else if (mpz_sgn(mpz) == -1 || mpz_sgn(rhs.mpz) == -1) {
      mpz_fdiv_q(mpz, mpz, rhs.mpz);
    } else if (mpz_sgn(mpz) == 0 || mpz_sgn(rhs.mpz) == 0) {
      mpz_set_ui(mpz, 0UL);
    } else {
      mpz_cdiv_q(mpz, mpz, rhs.mpz);
    }

    return *this;
  }

  inline BigInteger& operator/=(BIG_INTEGER_WORD rhs) {
    if (mpz_sgn(mpz) == -1) {
      (void) mpz_fdiv_q_ui(mpz, mpz, rhs);
    } else if (mpz_sgn(mpz) == 0 || rhs == 0) {
      mpz_set_ui(mpz, 0UL);
    } else {
      (void) mpz_cdiv_q_ui(mpz, mpz, rhs);
    }

    return *this;
  }

  inline BigInteger& operator/=(BIG_INTEGER_HALF_WORD rhs) {
    if (mpz_sgn(mpz) == -1) {
      (void) mpz_fdiv_q_ui(mpz, mpz, rhs);
    } else if (mpz_sgn(mpz) == 0 || rhs == 0) {
      mpz_set_ui(mpz, 0UL);
    } else {
      (void) mpz_cdiv_q_ui(mpz, mpz, rhs);
    }

    return *this;
  }

  inline BigInteger& operator%=(const BigInteger& rhs) {
    if (mpz_sgn(mpz) == -1 && mpz_sgn(rhs.mpz) == -1) {
      (void) mpz_cdiv_r(mpz, mpz, rhs.mpz);
    } else if (mpz_sgn(mpz) == -1 || mpz_sgn(rhs.mpz) == -1) {
      (void) mpz_fdiv_r(mpz, mpz, rhs.mpz);
    } else if (mpz_sgn(mpz) == 0 || mpz_sgn(rhs.mpz) == 0) {
      mpz_set_ui(mpz, 0UL);
    } else {
      (void) mpz_cdiv_r(mpz, mpz, rhs.mpz);
    }

    return *this;
  }

  inline BigInteger& operator%=(BIG_INTEGER_WORD rhs) {
    if (mpz_sgn(mpz) == -1) {
      (void) mpz_fdiv_r_ui(mpz, mpz, rhs);
    } else if (mpz_sgn(mpz) == 0 || rhs == 0) {
      mpz_set_ui(mpz, 0UL);
    } else {
      (void) mpz_cdiv_r_ui(mpz, mpz, rhs);
    }

    return *this;
  }

  inline BigInteger& operator%=(BIG_INTEGER_HALF_WORD rhs) {
    return BigInteger::operator%=(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  template<typename __ST>
  inline BigInteger& operator<<=(__ST sq) {
    static_assert(std::is_integral<__ST>::value,
                  "shift quantity must be an integral type!");
    static_assert(!std::is_signed<__ST>::value,
                  "shift quantity must be an unsigned integral type!");

    if (!sq)
      return *this;

    mpz_t tmp;
    mpz_init2(tmp, MPZ_INTEGER_BITS);

    for (__ST i = 0; i < sq + 1; ++i) {
      mpz_mul_2exp(tmp, mpz, i);
    }

    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    return *this;
  }

  template<typename __ST>
  inline BigInteger& operator>>=(__ST sq) {
    static_assert(std::is_integral<__ST>::value,
                  "shift quantity must be an integral type!");
    static_assert(!std::is_signed<__ST>::value,
                  "shift quantity must be an unsigned integral type!");

    if (!sq)
      return *this;

    mpz_t tmp;
    mpz_init2(tmp, MPZ_INTEGER_BITS);

    for (__ST i = 0; i < sq + 1; ++i) {
      mpz_mul_2exp(tmp, mpz, i);
      mpz_fdiv_q_2exp(tmp, mpz, i);
    }

    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    return *this;
  }

  inline BigInteger& operator&=(const BigInteger& rhs) {
    mpz_t tmp;
    mpz_init2(tmp, MPZ_INTEGER_BITS);
    mpz_and(tmp, mpz, rhs.mpz);
    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    return *this;
  }

  inline BigInteger& operator^=(const BigInteger& rhs) {
    mpz_t tmp;
    mpz_init2(tmp, MPZ_INTEGER_BITS);
    mpz_xor(tmp, mpz, rhs.mpz);
    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    return *this;
  }

  inline BigInteger& operator|=(const BigInteger& rhs) {
    mpz_t tmp;
    mpz_init2(tmp, MPZ_INTEGER_BITS);
    mpz_ior(tmp, mpz, rhs.mpz);
    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    return *this;
  }

  inline BigInteger& operator&=(BIG_INTEGER_WORD rhs) {
    mpz_t tmp;
    mpz_t smpz;
    mpz_init2(tmp, MPZ_INTEGER_BITS);
    mpz_init2(smpz, MPZ_INTEGER_BITS);
    mpz_set_ui(smpz, rhs);
    mpz_and(tmp, mpz, smpz);
    mpz_set(mpz, tmp);
    mpz_clear(smpz);
    mpz_clear(tmp);
    return *this;
  }

  inline BigInteger& operator^=(BIG_INTEGER_WORD rhs) {
    mpz_t tmp;
    mpz_t smpz;
    mpz_init2(tmp, MPZ_INTEGER_BITS);
    mpz_init2(smpz, MPZ_INTEGER_BITS);
    mpz_set_ui(smpz, rhs);
    mpz_xor(tmp, mpz, smpz);
    mpz_set(mpz, tmp);
    mpz_clear(smpz);
    mpz_clear(tmp);
    return *this;
  }

  inline BigInteger& operator|=(BIG_INTEGER_WORD rhs) {
    mpz_t tmp;
    mpz_t smpz;
    mpz_init2(tmp, MPZ_INTEGER_BITS);
    mpz_init2(smpz, MPZ_INTEGER_BITS);
    mpz_set_ui(smpz, rhs);
    mpz_ior(tmp, mpz, smpz);
    mpz_set(mpz, tmp);
    mpz_clear(smpz);
    mpz_clear(tmp);
    return *this;
  }

  inline void increment(BIG_INTEGER_WORD val) {
    mpz_add_ui(mpz, mpz, val);
  }

  inline void increment(BIG_INTEGER_WORD val) const {
    mpz_add_ui(mpz, mpz, val);
  }

  inline void decrement(BIG_INTEGER_WORD val) {
    mpz_sub_ui(mpz, mpz, val);
  }

  inline void decrement(BIG_INTEGER_WORD val) const {
    mpz_sub_ui(mpz, mpz, val);
  }

  inline bool is_zero() const {
    return mpz_sgn(mpz) == 0;
  }

  inline bool is_negative() const {
    return mpz_sgn(mpz) < 0;
  }

  inline bool is_positive() const {
    return mpz_sgn(mpz) > 0;
  }

  static inline int32_t log2(const BigInteger& rhs) {
    BigInteger bip = rhs >> 1UL;
    if (bip.is_negative() || bip.is_zero()) {
      return -1;
    }

    int32_t r = 0;

    while (bip > 0UL) {
      bip >>= 1UL;
      ++r;
    }

    return r;
  }

  uint64_t as_unsigned_long() const {
    return mpz_get_ui(mpz);
  }

  double as_double() const {
    return mpz_get_d(mpz);
  }

  void set_zero() {
    mpz_set_ui(mpz, 0UL);
  }

  void set_zero() const {
    mpz_set_ui(mpz, 0UL);
  }
};

inline void bi_set_0(BigInteger* p) {
  p->set_zero();
}

inline BigInteger bi_copy(const BigInteger& in) {
  return BigInteger(in);
}

inline void bi_copy_ip(const BigInteger& in, BigInteger* out) {
  *out = in;
}

inline int bi_compare(const BigInteger& left, const BigInteger& right) {
  return mpz_cmp(left.mpz, right.mpz);
}

inline int bi_compare(const BigInteger& left, BIG_INTEGER_WORD right) {
  return mpz_cmp_ui(left.mpz, right);
}

inline int bi_compare_0(const BigInteger& left) {
  return left.is_zero() ? 1 : 0;
}

inline int bi_compare_1(const BigInteger& left) {
  if (left.is_zero())
    return 0;

  return left.is_negative() ? -1 : 1;
}

inline void bi_add_ip(BigInteger* left, const BigInteger& right) {
  *left += right;
}

inline void bi_sub_ip(BigInteger* left, const BigInteger& right) {
  *left -= right;
}

inline void bi_increment(BigInteger* pBigInt, BIG_INTEGER_WORD value) {
  *pBigInt += value;
}

inline void bi_decrement(BigInteger* pBigInt, BIG_INTEGER_WORD value) {
  *pBigInt -= value;
}

inline BigInteger bi_load(BIG_INTEGER_WORD* a) {
  return BigInteger(*a);
}

inline BigInteger
bi_lshift_word(const BigInteger& left, BIG_INTEGER_WORD rightMult) {
  if (!rightMult)
    return left;

  return left << rightMult;
}

inline void bi_lshift_word_ip(BigInteger* left, BIG_INTEGER_WORD rightMult) {
  rightMult &= 63UL;
  if (!rightMult)
    return;

  BigInteger tmp = *left;
  *left = tmp << rightMult;
}

inline BigInteger
bi_rshift_word(const BigInteger& left, BIG_INTEGER_WORD rightMult) {
  if (!rightMult)
    return left;

  return left >> rightMult;
}

inline void bi_rshift_word_ip(BigInteger* left, BIG_INTEGER_WORD rightMult) {
  if (!rightMult)
    return;

  BigInteger tmp = *left;
  *left = tmp >> rightMult;
}

inline void bi_lshift_ip(BigInteger* left, BIG_INTEGER_WORD right) {
  BigInteger tmp = *left;
  *left = tmp << right;
}

inline void bi_rshift_ip(BigInteger* left, BIG_INTEGER_WORD right) {
  BigInteger tmp = *left >> right;
  *left = tmp;
}

inline int bi_log2(const BigInteger& n) {
  int pw = 0;
  BigInteger p = n >> 1UL;
  while (bi_compare_0(p) != 0) {
    bi_rshift_ip(&p, 1U);
    ++pw;
  }

  return pw;
}

inline int bi_and_1(const BigInteger& left) {
  return mpz_tstbit(left.mpz, 0);
}

inline void bi_and_ip(BigInteger* left, const BigInteger& right) {
  *left &= right;
}

inline void bi_or_ip(BigInteger* left, const BigInteger& right) {
  *left |= right;
}

inline void bi_xor_ip(BigInteger* left, const BigInteger& right) {
  *left ^= right;
}

inline void bi_not_ip(BigInteger* left) {
  *left = ~(*left);
}

inline double bi_to_double(const BigInteger& in) {
  return in.operator double();
}

void bi_div_mod_small(const BigInteger& left, BIG_INTEGER_HALF_WORD right,
                      BigInteger* quotient, BIG_INTEGER_HALF_WORD* rmndr);

void bi_div_mod(const BigInteger& left, const BigInteger& right,
                BigInteger* quotient, BigInteger* rmndr);

std::ostream& operator<<(std::ostream& os, const BigInteger& bi);

} // namespace Qrack

