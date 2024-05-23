//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qimcifa contributors, 2022, 2023. All rights reserved.
// (C) Quantum Circuits, Inc, 2024. All rights reserved.
//
// This header has been adapted for OpenCL and C, from big_integer.c by Andre Azevedo.
// This file has been restored to sanity by Stefan Teleman for Quantum Circuits, Inc.
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

#include <stdexcept>

#include "big_integer.hpp"

static void (*gmp_free_mem_func)(void*, size_t);

namespace Qrack {

BigInteger operator*(const BigInteger& left, BIG_INTEGER_HALF_WORD right) {
  BigInteger result;
  mpz_mul_ui(result.mpz, left.mpz, right);
  return result;
}

void bi_div_mod_small(const BigInteger& left, BIG_INTEGER_HALF_WORD right,
                      BigInteger* quotient, BIG_INTEGER_HALF_WORD* rmndr) {
  if (left.is_zero()) {
    mpz_set_ui(quotient->mpz, 0UL);
    *rmndr = 0U;
  } else if (left.is_negative()) {
    mpz_t r;
    mpz_init2(r, BIG_INTEGER_WORD_BITS);
    *rmndr = static_cast<BIG_INTEGER_HALF_WORD>(mpz_fdiv_qr_ui(quotient->mpz, r,
                                                               left.mpz, right));
    mpz_clear(r);
  } else {
    mpz_t r;
    mpz_init2(r, BIG_INTEGER_WORD_BITS);
    *rmndr = static_cast<BIG_INTEGER_HALF_WORD>(mpz_cdiv_qr_ui(quotient->mpz, r,
                                                               left.mpz, right));
    mpz_clear(r);
  }
}

void bi_div_mod(const BigInteger& left, const BigInteger& right,
                BigInteger* quotient, BigInteger* rmndr) {
  if (quotient && rmndr) {
    mpz_cdiv_qr(quotient->mpz, rmndr->mpz, left.mpz, right.mpz);
  } else if (quotient && !rmndr) {
    mpz_cdiv_q(quotient->mpz, left.mpz, right.mpz);
  } else if (rmndr && !quotient) {
    mpz_cdiv_r(rmndr->mpz, left.mpz, right.mpz);
  }
}

std::ostream& operator<<(std::ostream& os, const BigInteger& bi) {
  if (char* bip = mpz_get_str(NULL, 10, bi.mpz)) {
    std::string bis = bip;
    os << bis;
    mp_get_memory_functions(NULL, NULL, &gmp_free_mem_func);
    gmp_free_mem_func(bip, std::strlen(bip) + 1);
  }

  return os;
}

} // namespace Qrack

