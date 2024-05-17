
//
// (C) Daniel Strano and the Qimcifa contributors, 2022, 2023. All rights reserved.
//
// This header has been adapted for OpenCL and C, from big_integer.c by Andre Azevedo.
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

#include <gmp.h>

#define BIG_INTEGER_WORD_BITS 64U
#define BIG_INTEGER_WORD_POWER 6U
#define BIG_INTEGER_WORD uint64_t
#define BIG_INTEGER_HALF_WORD uint32_t
#define BIG_INTEGER_HALF_WORD_POW 0x100000000ULL
#define BIG_INTEGER_HALF_WORD_MASK 0xFFFFFFFFULL
#define BIG_INTEGER_HALF_WORD_MASK_NOT 0xFFFFFFFF00000000ULL

// This can be any power of 2 greater than (or equal to) 64:
#define BIG_INTEGER_BITS (1 << QBCAPPOW)
#define BIG_INTEGER_WORD_SIZE (long long)(BIG_INTEGER_BITS / BIG_INTEGER_WORD_BITS)

// The rest of the constants need to be consistent with the one above:
constexpr size_t BIG_INTEGER_HALF_WORD_BITS = BIG_INTEGER_WORD_BITS >> 1U;
constexpr int BIG_INTEGER_HALF_WORD_SIZE = BIG_INTEGER_WORD_SIZE << 1U;
constexpr int BIG_INTEGER_MAX_WORD_INDEX = BIG_INTEGER_WORD_SIZE - 1U;

class BigInteger {
public:
  inline void sync_bits() {
    size_t sz = mpz_sizeinbase(mpz, 2) + CHAR_BIT - 1;
    void* p = mpz_export(reinterpret_cast<void*>(bits), &sz, 1, 1, 0, 0, mpz);
    (void) p;
    assert(p == bits && "invalid pointer obtained from mpz-export!");
  }

  inline void sync_bits() const {
    size_t sz = mpz_sizeinbase(mpz, 2) + CHAR_BIT - 1;
    void* p = mpz_export(reinterpret_cast<void*>(bits), &sz, 1, 1, 0, 0, mpz);
    (void) p;
    assert(p == bits && "invalid pointer obtained from mpz-export!");
  }

  inline void sync_bits(const mpz_t& impz) {
    size_t sz = mpz_sizeinbase(impz, 2) + CHAR_BIT - 1;
    void* p = mpz_export(reinterpret_cast<void*>(bits), &sz, 1, 1, 0, 0, impz);
    (void) p;
    assert(p == bits && "invalid pointer obtained from mpz-export!");
  }

  inline void sync_bits(const mpz_t& impz) const {
    size_t sz = mpz_sizeinbase(impz, 2) + CHAR_BIT - 1;
    void* p = mpz_export(reinterpret_cast<void*>(bits), &sz, 1, 1, 0, 0, impz);
    (void) p;
    assert(p == bits && "invalid pointer obtained from mpz-export!");
  }

  inline void sync_mpz() {
    mpz_import(mpz, sizeof(bits), 1, 1, 0, 0, &bits[0]);
  }

  inline void sync_mpz() const {
    mpz_import(mpz, sizeof(bits), 1, 1, 0, 0, &bits[0]);
  }

  inline void sync_mpz(mpz_t& impz) {
    mpz_import(impz, sizeof(bits), 1, 1, 0, 0, &bits[0]);
  }

  inline void sync_mpz(mpz_t& impz) const {
    mpz_import(impz, sizeof(bits), 1, 1, 0, 0, &bits[0]);
  }

public:
  mutable BIG_INTEGER_WORD bits[BIG_INTEGER_WORD_SIZE];
  mutable mpz_t mpz;

public:
  BigInteger() : bits{0}, mpz() {
    mpz_init2(mpz, BIG_INTEGER_WORD_BITS);
  }

  BigInteger(const BigInteger& rhs) : bits{0}, mpz() {
    mpz_init2(mpz, BIG_INTEGER_WORD_BITS);
    mpz_set(mpz, rhs.mpz);
    sync_bits(rhs.mpz);
  }

  BigInteger(BigInteger&& rhs) : bits{0}, mpz() {
    mpz_init2(mpz, BIG_INTEGER_WORD_BITS);
    mpz_set(mpz, rhs.mpz);
    sync_bits(rhs.mpz);
    mpz_set_ui(rhs.mpz, 0UL);
  }

  BigInteger(BIG_INTEGER_WORD rhs) : bits{0}, mpz() {
    bits[0] = rhs;
    mpz_init2(mpz, BIG_INTEGER_WORD_BITS);
    sync_mpz();
  }

  BigInteger(const mpz_t& rhs) : bits{0}, mpz() {
    mpz_init2(mpz, BIG_INTEGER_WORD_BITS);
    mpz_set(mpz, rhs);
    sync_bits(rhs);
  }

  BigInteger(mpz_t&& rhs) : bits{0}, mpz() {
    mpz_init2(mpz, BIG_INTEGER_WORD_BITS);
    mpz_set(mpz, rhs);
    sync_bits(rhs);
    mpz_set_ui(rhs, 0UL);
  }

  ~BigInteger() {
    (void) mpz_clear(mpz);
  }

  BigInteger& operator=(const BigInteger& rhs) {
    if (this != &rhs) {
      mpz_set(mpz, rhs.mpz);
      sync_bits();
    }

    return *this;
  }

  BigInteger& operator=(BigInteger&& rhs) {
    if (this != &rhs) {
      mpz_set(mpz, rhs.mpz);
      sync_bits();
    }

    return *this;
  }

  BigInteger& operator=(BIG_INTEGER_WORD rhs) {
    mpz_set_ui(mpz, rhs);
    sync_bits();
    return *this;
  }

  BigInteger& operator=(const mpz_t& rhs) {
    mpz_set(mpz, rhs);
    sync_bits();
    return *this;
  }

  explicit inline operator BIG_INTEGER_WORD() {
    sync_bits();
    return bits[0U];
  }

  explicit inline operator BIG_INTEGER_WORD() const {
    sync_bits();
    return bits[0U];
  }

  explicit inline operator uint32_t() {
    sync_bits();
    uint32_t* p32 = reinterpret_cast<uint32_t*>(bits);
    return p32[0U];
  }

  explicit inline operator uint32_t() const {
    sync_bits();
    const uint32_t* p32 = reinterpret_cast<const uint32_t*>(bits);
    return p32[0U];
  }

  explicit inline operator uint16_t() {
    sync_bits();
    uint16_t* p16 = reinterpret_cast<uint16_t*>(bits);
    return p16[0U];
  }

  explicit inline operator uint16_t() const {
    sync_bits();
    const uint16_t* p16 = reinterpret_cast<const uint16_t*>(bits);
    return p16[0U];
  }

  explicit inline operator uint8_t() {
    sync_bits();
    uint8_t* p8 = reinterpret_cast<uint8_t*>(bits);
    return p8[0U];
  }

  explicit inline operator uint8_t() const {
    sync_bits();
    const uint8_t* p8 = reinterpret_cast<const uint8_t*>(bits);
    return p8[0U];
  }

  // Standard Arithmetic Operators.
  inline BigInteger operator+(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_add(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator+(BIG_INTEGER_WORD rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_add_ui(result, mpz, rhs);
    return result;
  }

  inline BigInteger operator+(BIG_INTEGER_HALF_WORD rhs) const {
    return BigInteger::operator+(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger operator-(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_sub(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator-(BIG_INTEGER_WORD rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_sub_ui(result, mpz, rhs);
    return result;
  }

  inline BigInteger operator-(BIG_INTEGER_HALF_WORD rhs) const {
    return BigInteger::operator-(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger operator*(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_mul(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator*(BIG_INTEGER_WORD rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_mul_ui(result, mpz, rhs);
    return result;
  }

  inline BigInteger operator*(BIG_INTEGER_HALF_WORD rhs) const {
    return BigInteger::operator*(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger operator/(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);

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
    mpz_init2(result, BIG_INTEGER_WORD_BITS);

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
    mpz_init2(result, BIG_INTEGER_WORD_BITS);

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
    mpz_init2(result, BIG_INTEGER_WORD_BITS);

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
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
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

    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);

    for (uint32_t i = 0; i < sq + 1U; ++i) {
      mpz_mul_2exp(result, mpz, i);
      mpz_tdiv_r_2exp(result, result, i);
    }

    return result;
  }

  template<typename __ST>
  inline BigInteger operator>>(__ST sq) const {
    static_assert(std::is_integral<__ST>::value,
                  "shift quantity must be an integral type!");
    static_assert(!std::is_signed<__ST>::value,
                  "shift quantity must be an unsigned integral type!");

    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);

    for (uint32_t i = 0; i < sq + 1U; ++i) {
      mpz_mul_2exp(result, mpz, i);
      mpz_fdiv_q_2exp(result, mpz, i);
    }

    return result;
  }

  inline BigInteger operator&(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_and(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator^(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_xor(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator|(const BigInteger& rhs) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_ior(result, mpz, rhs.mpz);
    return result;
  }

  inline BigInteger operator&(BIG_INTEGER_WORD rhs) const {
    mpz_t op;
    mpz_t result;
    mpz_init2(op, BIG_INTEGER_WORD_BITS);
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_set_ui(op, rhs);
    mpz_and(result, mpz, op);
    mpz_clear(op);
    return result;
  }

  inline BigInteger operator^(BIG_INTEGER_WORD rhs) const {
    mpz_t op;
    mpz_t result;
    mpz_init2(op, BIG_INTEGER_WORD_BITS);
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_set_ui(op, rhs);
    mpz_xor(result, mpz, op);
    mpz_clear(op);
    return result;
  }

  inline BigInteger operator|(BIG_INTEGER_WORD rhs) const {
    mpz_t op;
    mpz_t result;
    mpz_init2(op, BIG_INTEGER_WORD_BITS);
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_set_ui(op, rhs);
    mpz_ior(result, mpz, op);
    mpz_clear(op);
    return result;
  }

  inline BigInteger& operator++() {
    mpz_add_ui(mpz, mpz, 1UL);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator--() {
    mpz_sub_ui(mpz, mpz, 1UL);
    sync_bits();
    return *this;
  }

  inline BigInteger operator++(int) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_set(result, mpz);
    mpz_add_ui(mpz, mpz, 1UL);
    return result;
  }

  inline BigInteger operator--(int) const {
    mpz_t result;
    mpz_init2(result, BIG_INTEGER_WORD_BITS);
    mpz_set(result, mpz);
    mpz_sub_ui(mpz, mpz, 1UL);
    return result;
  }

  inline BigInteger& operator+=(const BigInteger& rhs) {
    mpz_add(mpz, rhs.mpz, mpz);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator+=(BIG_INTEGER_WORD rhs) {
    mpz_add_ui(mpz, mpz, rhs);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator+=(BIG_INTEGER_HALF_WORD rhs) {
    return BigInteger::operator+=(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger& operator-=(const BigInteger& rhs) {
    mpz_sub(mpz, mpz, rhs.mpz);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator-=(BIG_INTEGER_WORD rhs) {
    mpz_sub_ui(mpz, mpz, rhs);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator-=(BIG_INTEGER_HALF_WORD rhs) {
    return BigInteger::operator-=(static_cast<BIG_INTEGER_WORD>(rhs));
  }

  inline BigInteger& operator*=(const BigInteger& rhs) {
    mpz_mul(mpz, mpz, rhs.mpz);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator*=(BIG_INTEGER_WORD rhs) {
    mpz_mul_ui(mpz, mpz, rhs);
    sync_bits();
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

    sync_bits();
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

    sync_bits();
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

    sync_bits();
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

    sync_bits();
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

    sync_bits();
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

    mpz_t tmp;
    mpz_init2(tmp, BIG_INTEGER_WORD_BITS * 2U);

    for (uint32_t i = 0; i < sq + 1; ++i) {
      mpz_mul_2exp(tmp, mpz, i);
      mpz_tdiv_r_2exp(tmp, tmp, i);
    }

    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    sync_bits();
    return *this;
  }

  template<typename __ST>
  inline BigInteger& operator>>=(__ST sq) {
    static_assert(std::is_integral<__ST>::value,
                  "shift quantity must be an integral type!");
    static_assert(!std::is_signed<__ST>::value,
                  "shift quantity must be an unsigned integral type!");

    mpz_t tmp;
    mpz_init2(tmp, BIG_INTEGER_WORD_BITS);

    for (uint32_t i = 0; i < sq + 1; ++i) {
      mpz_mul_2exp(tmp, mpz, i);
      mpz_fdiv_q_2exp(tmp, mpz, i);
    }

    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator&=(const BigInteger& rhs) {
    mpz_t tmp;
    mpz_init2(tmp, BIG_INTEGER_WORD_BITS);
    mpz_and(tmp, mpz, rhs.mpz);
    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator^=(const BigInteger& rhs) {
    mpz_t tmp;
    mpz_init2(tmp, BIG_INTEGER_WORD_BITS);
    mpz_xor(tmp, mpz, rhs.mpz);
    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator|=(const BigInteger& rhs) {
    mpz_t tmp;
    mpz_init2(tmp, BIG_INTEGER_WORD_BITS);
    mpz_ior(tmp, mpz, rhs.mpz);
    mpz_set(mpz, tmp);
    mpz_clear(tmp);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator&=(BIG_INTEGER_WORD rhs) {
    mpz_t tmp;
    mpz_t smpz;
    mpz_init2(tmp, BIG_INTEGER_WORD_BITS);
    mpz_init2(smpz, BIG_INTEGER_WORD_BITS);
    mpz_set_ui(smpz, rhs);
    mpz_and(tmp, mpz, smpz);
    mpz_set(mpz, tmp);
    mpz_clear(smpz);
    mpz_clear(tmp);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator^=(BIG_INTEGER_WORD rhs) {
    mpz_t tmp;
    mpz_t smpz;
    mpz_init2(tmp, BIG_INTEGER_WORD_BITS);
    mpz_init2(smpz, BIG_INTEGER_WORD_BITS);
    mpz_set_ui(smpz, rhs);
    mpz_xor(tmp, mpz, smpz);
    mpz_set(mpz, tmp);
    mpz_clear(smpz);
    mpz_clear(tmp);
    sync_bits();
    return *this;
  }

  inline BigInteger& operator|=(BIG_INTEGER_WORD rhs) {
    mpz_t tmp;
    mpz_t smpz;
    mpz_init2(tmp, BIG_INTEGER_WORD_BITS);
    mpz_init2(smpz, BIG_INTEGER_WORD_BITS);
    mpz_set_ui(smpz, rhs);
    mpz_ior(tmp, mpz, smpz);
    mpz_set(mpz, tmp);
    mpz_clear(smpz);
    mpz_clear(tmp);
    sync_bits();
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
    BigInteger bip = rhs >> 1U;
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
};

inline void bi_set_0(BigInteger* p) {
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    p->bits[i] = 0U;
  }

  p->sync_mpz();
}

inline BigInteger bi_copy(const BigInteger& in) {
  BigInteger result;
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    result.bits[i] = in.bits[i];
  }

  result.sync_mpz();
  return result;
}

inline void bi_copy_ip(const BigInteger& in, BigInteger* out) {
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    out->bits[i] = in.bits[i];
  }
}

inline int bi_compare(const BigInteger& left, const BigInteger& right) {
  return mpz_cmp(left.mpz, right.mpz);

#if 0
  // ORIGINAL:
  for (int i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
    if (left.bits[i] > right.bits[i]) {
      return 1;
    }
    if (left.bits[i] < right.bits[i]) {
      return -1;
    }
  }

  return 0;
#endif
}

inline int bi_compare_0(const BigInteger& left) {
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    if (left.bits[i]) {
      return 1;
    }
  }

  return 0;
}

inline int bi_compare_1(const BigInteger& left) {
  for (int i = BIG_INTEGER_MAX_WORD_INDEX; i > 0; --i) {
    if (left.bits[i]) {
      return 1;
    }
  }

  if (left.bits[0] > 1U) {
    return 1;
  }
  if (left.bits[0] < 1U) {
    return -1;
  }

  return 0;
}

#if 0
inline BigInteger operator+(const BigInteger& left, const BigInteger& right)
{
    BigInteger result;
    result.bits[0U] = 0U;
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result.bits[i] += left.bits[i] + right.bits[i];
        result.bits[i + 1] = (result.bits[i] < left.bits[i]) ? 1 : 0;
    }
    result.bits[BIG_INTEGER_MAX_WORD_INDEX] += right.bits[BIG_INTEGER_MAX_WORD_INDEX];

    return result;
}
#endif

inline void bi_add_ip(BigInteger* left, const BigInteger& right) {
  for (uint32_t i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
    BIG_INTEGER_WORD temp = left->bits[i];
    left->bits[i] += right.bits[i];
    uint32_t j = i;
    while ((j < BIG_INTEGER_MAX_WORD_INDEX) && (left->bits[j] < temp)) {
      temp = left->bits[++j]++;
    }
  }

  left->bits[BIG_INTEGER_MAX_WORD_INDEX] += right.bits[BIG_INTEGER_MAX_WORD_INDEX];
  left->sync_mpz();
}

#if 0
inline BigInteger operator-(const BigInteger& left, const BigInteger& right)
{
    BigInteger result;
    result.bits[0U] = 0U;
    for (int i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
        result.bits[i] += left.bits[i] - right.bits[i];
        result.bits[i + 1] = (result.bits[i] > left.bits[i]) ? -1 : 0;
    }
    result.bits[BIG_INTEGER_MAX_WORD_INDEX] -= right.bits[BIG_INTEGER_MAX_WORD_INDEX];

    return result;
}
#endif

inline void bi_sub_ip(BigInteger* left, const BigInteger& right) {
  for (uint32_t i = 0; i < BIG_INTEGER_MAX_WORD_INDEX; ++i) {
    BIG_INTEGER_WORD temp = left->bits[i];
    left->bits[i] -= right.bits[i];
    uint32_t j = i;
    while ((j < BIG_INTEGER_MAX_WORD_INDEX) && (left->bits[j] > temp)) {
      temp = left->bits[++j]--;
    }
  }

  left->bits[BIG_INTEGER_MAX_WORD_INDEX] -= right.bits[BIG_INTEGER_MAX_WORD_INDEX];
  left->sync_mpz();
}

inline void bi_increment(BigInteger* pBigInt, const BIG_INTEGER_WORD& value) {
  BIG_INTEGER_WORD temp = pBigInt->bits[0];
  pBigInt->bits[0] += value;
  pBigInt->sync_mpz();

  if (temp <= pBigInt->bits[0])
    return;

  for (uint32_t i = 1; i < BIG_INTEGER_WORD_SIZE; i++) {
    temp = pBigInt->bits[i]++;
    if (temp <= pBigInt->bits[i]) {
      break;
    }
  }

  pBigInt->sync_mpz();
}

inline void bi_decrement(BigInteger* pBigInt, const BIG_INTEGER_WORD& value) {
  BIG_INTEGER_WORD temp = pBigInt->bits[0];
  pBigInt->bits[0] -= value;
  pBigInt->sync_mpz();

  if (temp >= pBigInt->bits[0])
    return;

  for (uint32_t i = 0; i < BIG_INTEGER_WORD_SIZE; i++) {
    temp = pBigInt->bits[i]--;
    if (temp >= pBigInt->bits[i]) {
      break;
    }
  }

  pBigInt->sync_mpz();
}

inline BigInteger bi_load(BIG_INTEGER_WORD* a) {
  BigInteger result;
  for (uint32_t i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    result.bits[i] = a[i];
  }

  result.sync_mpz();
  return result;
}

inline BigInteger
bi_lshift_word(const BigInteger& left, BIG_INTEGER_WORD rightMult) {
  if (!rightMult) {
    return left;
  }

  BigInteger result = 0UL;
  for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
    result.bits[i] = left.bits[i - rightMult];
  }

  result.sync_mpz();
  return result;
}

inline void bi_lshift_word_ip(BigInteger* left, BIG_INTEGER_WORD rightMult) {
  rightMult &= 63U;
  if (!rightMult)
    return;

  for (uint32_t i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
    left->bits[i] = left->bits[i - rightMult];
  }

  for (BIG_INTEGER_WORD i = 0UL; i < rightMult; ++i) {
    left->bits[i] = 0U;
  }

  left->sync_mpz();
}

inline BigInteger
bi_rshift_word(const BigInteger& left, const BIG_INTEGER_WORD& rightMult) {
  if (!rightMult)
    return left;

  BigInteger result = 0UL;
  for (int i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
    result.bits[i - rightMult] = left.bits[i];
  }

  result.sync_mpz();
  return result;
}

inline void bi_rshift_word_ip(BigInteger* left, const BIG_INTEGER_WORD& rightMult) {
  if (!rightMult)
    return;

  for (uint32_t i = rightMult; i < BIG_INTEGER_WORD_SIZE; ++i) {
    left->bits[i - rightMult] = left->bits[i];
  }

  for (BIG_INTEGER_WORD i = 0UL; i < rightMult; ++i) {
    left->bits[BIG_INTEGER_MAX_WORD_INDEX - i] = 0UL;
  }

  left->sync_mpz();
}

inline BigInteger operator<<(const BigInteger& left, BIG_INTEGER_WORD right) {
  const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
  const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

  BigInteger result = bi_lshift_word(left, rShift64);
  if (!rMod)
    return result;

  const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
  BIG_INTEGER_WORD carry = 0U;
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    right = result.bits[i];
    result.bits[i] = carry | (right << rMod);
    carry = right >> rModComp;
  }

  result.sync_mpz();
  return result;
}

inline void bi_lshift_ip(BigInteger* left, BIG_INTEGER_WORD right) {
  const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
  const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

  bi_lshift_word_ip(left, rShift64);
  if (!rMod) {
    return;
  }

  const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
  BIG_INTEGER_WORD carry = 0U;
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    right = left->bits[i];
    left->bits[i] = carry | (right << rMod);
    carry = right >> rModComp;
  }

  left->sync_mpz();
}

inline BigInteger operator>>(const BigInteger& left, BIG_INTEGER_WORD right) {
  const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
  const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

  BigInteger result = bi_rshift_word(left, rShift64);
  if (!rMod) {
    return result;
  }

  const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
  BIG_INTEGER_WORD carry = 0U;

  for (int32_t i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
    right = result.bits[i];
    result.bits[i] = carry | (right >> rMod);
    carry = right << rModComp;
  }

  result.sync_mpz();
  return result;
}

inline void bi_rshift_ip(BigInteger* left, BIG_INTEGER_WORD right) {
  const int rShift64 = right >> BIG_INTEGER_WORD_POWER;
  const int rMod = right - (rShift64 << BIG_INTEGER_WORD_POWER);

  bi_rshift_word_ip(left, rShift64);
  if (!rMod)
    return;

  const int rModComp = BIG_INTEGER_WORD_BITS - rMod;
  BIG_INTEGER_WORD carry = 0U;

  for (int32_t i = BIG_INTEGER_MAX_WORD_INDEX; i >= 0; --i) {
    right = left->bits[i];
    left->bits[i] = carry | (right >> rMod);
    carry = right << rModComp;
  }

  left->sync_mpz();
}

inline int bi_log2(const BigInteger& n) {
  int pw = 0;
  BigInteger p = n >> 1U;
  while (bi_compare_0(p) != 0) {
    bi_rshift_ip(&p, 1U);
    ++pw;
  }

  return pw;
}

inline int bi_and_1(const BigInteger& left) { return left.bits[0] & 1; }

#if 0
inline BigInteger operator&(const BigInteger& left, const BigInteger& right) {
  return left & right;

#if 0
  // ORIGINAL:
  BigInteger result;
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    result.bits[i] = left.bits[i] & right.bits[i];
  }

  result.sync_mpz();
  return result;
#endif
}
#endif

inline void bi_and_ip(BigInteger* left, const BigInteger& right) {
  *left &= right;

#if 0
  // ORIGINAL:
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    left->bits[i] &= right.bits[i];
  }

  left->sync_mpz();
#endif
}

#if 0
inline BigInteger operator|(const BigInteger& left, const BigInteger& right) {
  return left | right;

#if 0
  // ORIGINAL:
  BigInteger result;
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    result.bits[i] = left.bits[i] | right.bits[i];
  }

  return result;
#endif
}
#endif

inline void bi_or_ip(BigInteger* left, const BigInteger& right) {
  *left |= right;

#if 0
  // ORIGINAL;
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    left->bits[i] |= right.bits[i];
  }
#endif
}

#if 0
inline BigInteger operator^(const BigInteger& left, const BigInteger& right) {
  return left ^ right;

#if 0
  // ORIGINAL:
  BigInteger result;
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    result.bits[i] = left.bits[i] ^ right.bits[i];
  }

  return result;
#endif
}
#endif

inline void bi_xor_ip(BigInteger* left, const BigInteger& right) {
  *left ^= right;

#if 0
  // ORIGINAL:
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    left->bits[i] ^= right.bits[i];
  }
#endif
}

#if 0
inline BigInteger operator~(const BigInteger& left)
{
    BigInteger result;
    for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
        result.bits[i] = ~(left.bits[i]);
    }

    return result;
}
#endif

inline void bi_not_ip(BigInteger* left) {
  *left = ~(*left);

#if 0
  // ORIGINAL
  for (int i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    left->bits[i] = ~(left->bits[i]);
  }
#endif
}

inline double bi_to_double(const BigInteger& in) {
  double toRet = 0.0;
  for (uint32_t i = 0; i < BIG_INTEGER_WORD_SIZE; ++i) {
    if (in.bits[i]) {
      toRet += in.bits[i] * std::pow(2.0, BIG_INTEGER_WORD_BITS * i);
    }
  }

  return toRet;
}

#if 0
inline bool operator<(const BigInteger& left, const BigInteger& right) { return bi_compare(left, right) < 0; }

/**
 * "Schoolbook multiplication" (on half words)
 * Complexity - O(x^2)
 */
BigInteger operator*(const BigInteger& left, BIG_INTEGER_HALF_WORD right);

#if BIG_INTEGER_BITS > 80
/**
 * Adapted from Qrack! (The fundamental algorithm was discovered before.)
 * Complexity - O(log)
 */
BigInteger operator*(const BigInteger& left, const BigInteger& right);
#else
/**
 * "Schoolbook multiplication" (on half words)
 * Complexity - O(x^2)
 */
BigInteger operator*(const BigInteger& left, const BigInteger& right);
#endif
#endif // if 0

/**
 * "Schoolbook division" (on half words)
 * Complexity - O(x^2)
 */
void bi_div_mod_small(const BigInteger& left, BIG_INTEGER_HALF_WORD right,
                      BigInteger* quotient, BIG_INTEGER_HALF_WORD* rmndr);

/**
 * Adapted from Qrack! (The fundamental algorithm was discovered before.)
 * Complexity - O(log)
 */
void bi_div_mod(const BigInteger& left, const BigInteger& right,
                BigInteger* quotient, BigInteger* rmndr);

