//////////////////////////////////////////////////////////////////////////////////////
//
// (C) Daniel Strano and the Qrack contributors 2017-2023. All rights reserved.
//
// This is a multithreaded, universal quantum register simulation, allowing
// (nonphysical) register cloning and direct measurement of probability and
// phase, to leverage what advantages classical emulation of qubits can have.
//
// Licensed under the GNU Lesser General Public License V3.
// See LICENSE.md in the project root or https://www.gnu.org/licenses/lgpl-3.0.en.html
// for details.

#pragma once

#define _USE_MATH_DEFINES

#include "config.h"

#include <complex>
#include <functional>
#include <climits>
#include <cmath>
#include <memory>
#include <cstring>
#include <cstdint>

#define IS_NORM_0(c) (norm(c) <= FP_NORM_EPSILON)
#define IS_SAME(c1, c2) (IS_NORM_0((c1) - (c2)))
#define IS_OPPOSITE(c1, c2) (IS_NORM_0((c1) + (c2)))

#if ENABLE_CUDA
#include <cuda_runtime.h>
#endif

#if (FPPOW < 5) && !defined(__arm__)
#include "half.hpp"
#endif

#if QBCAPPOW < 32U
using bitLenInt = uint32_t;
#else
using bitLenInt = uint64_t;
#endif

#if UINTPOW < 4
using bitCapIntOcl = uint8_t;
#elif UINTPOW < 5
using bitCapIntOcl = uint16_t;
#elif UINTPOW < 6
using bitCapIntOcl = uint32_t
#else
using bitCapIntOcl = uint64_t;
#endif

#if QBCAPPOW < 6
using bitCapInt = uint32_t;
#elif QBCAPPOW < 7
using bitCapInt = uint64_t;
#else
#include "big_integer.hpp"
using bitCapInt = Qrack::BigInteger;
#endif

#if FPPOW < 5
#ifdef __arm__
namespace Qrack {
using complex = std::complex<__fp16>;
typedef __fp16 real1;
typedef float real1_f;
typedef float real1_s;
#else
#if (CPP_STD >= 23) && __has_include(<stdfloat>)
#include <stdfloat>
#endif
#if defined(__STDCPP_FLOAT16_T__)
namespace Qrack {
using complex = std::complex<float16_t>;
typedef float16_t real1;
typedef float real1_f;
typedef float real1_s;
#else
#include "half.hpp"
namespace Qrack {
using complex = std::complex<half_float::half>;
typedef half_float::half real1;
typedef float real1_f;
typedef float real1_s;
#endif
#endif
#elif FPPOW < 6
namespace Qrack {
using complex = std::complex<float>;
typedef float real1;
typedef float real1_f;
typedef float real1_s;
#elif FPPOW < 7
namespace Qrack {
using complex = std::complex<double>;
typedef double real1;
typedef double real1_f;
typedef double real1_s;
#else
#if (CPP_STD >= 23) && __has_include(<stdfloat>)
#include <stdfloat>
#endif
#if defined(__STDCPP_FLOAT128_T__)
namespace Qrack {
using complex = std::complex<float128_t>;
typedef float128_t real1;
typedef float128_t real1_f;
typedef double real1_s;
#else
#include <boost/multiprecision/float128.hpp>
#include <quadmath.h>
namespace Qrack {
using complex = std::complex<boost::multiprecision::float128>;
typedef boost::multiprecision::float128 real1;
typedef boost::multiprecision::float128 real1_f;
typedef double real1_s;
#endif
#endif

const bitCapInt ONE_BCI = bitCapInt(1UL);
const bitCapInt ZERO_BCI = bitCapInt(0UL);
constexpr bitLenInt bitsInCap = bitLenInt(1U) << bitLenInt(QBCAPPOW);

typedef std::shared_ptr<complex> BitOp;

// Called once per value between begin and end.
typedef std::function<void(const bitCapIntOcl&, const unsigned& cpu)> ParallelFunc;
typedef std::function<bitCapIntOcl(const bitCapIntOcl&)> IncrementFunc;
typedef std::function<bitCapInt(const bitCapInt&)> BdtFunc;
typedef std::function<void(const bitCapInt&, const unsigned& cpu)> ParallelFuncBdt;

class StateVector;
class StateVectorArray;
class StateVectorSparse;

typedef std::shared_ptr<StateVector> StateVectorPtr;
typedef std::shared_ptr<StateVectorArray> StateVectorArrayPtr;
typedef std::shared_ptr<StateVectorSparse> StateVectorSparsePtr;

typedef std::function<void(void)> DispatchFn;

class QEngine;
typedef std::shared_ptr<QEngine> QEnginePtr;

#define bitsInByte 8U
#define qrack_rand_gen std::mt19937_64
#define qrack_rand_gen_ptr std::shared_ptr<qrack_rand_gen>
#define QRACK_ALIGN_SIZE 64U

#if FPPOW < 5
#define QRACK_CONST const
const real1 ZERO_R1 = (real1)0.0f;
constexpr real1_f ZERO_R1_F = 0.0f;
const real1 ONE_R1 = (real1)1.0f;
constexpr real1_f ONE_R1_F = 1.0f;
const real1 REAL1_DEFAULT_ARG = (real1)-999.0f;
// Half the probability in any single permutation of 20 maximally superposed qubits
const real1 REAL1_EPSILON = (real1)0.000000477f;
const real1 PI_R1 = (real1)M_PI;
const real1 SQRT2_R1 = (real1)M_SQRT2;
const real1 SQRT1_2_R1 = (real1)M_SQRT1_2;
#elif FPPOW < 6
#define QRACK_CONST constexpr
#define ZERO_R1 0.0f
#define ZERO_R1_F 0.0f
#define ONE_R1 1.0f
#define ONE_R1_F 1.0f
constexpr real1 PI_R1 = (real1)M_PI;
constexpr real1 SQRT2_R1 = (real1)M_SQRT2;
constexpr real1 SQRT1_2_R1 = (real1)M_SQRT1_2;
#define REAL1_DEFAULT_ARG -999.0f
// Half the probability in any single permutation of 48 maximally superposed qubits
#define REAL1_EPSILON 1.7763568394002505e-15f
#elif FPPOW < 7
#define QRACK_CONST constexpr
#define ZERO_R1 0.0
#define ZERO_R1_F 0.0
#define ONE_R1 1.0
#define ONE_R1_F 1.0
#define PI_R1 M_PI
#define SQRT2_R1 M_SQRT2
#define SQRT1_2_R1 M_SQRT1_2
#define REAL1_DEFAULT_ARG -999.0
// Half the probability in any single permutation of 96 maximally superposed qubits
#define REAL1_EPSILON 6.310887241768095e-30
#else
#define QRACK_CONST constexpr
constexpr real1 ZERO_R1 = (real1)0.0;
#define ZERO_R1_F 0.0
constexpr real1 ONE_R1 = (real1)1.0;
#define ONE_R1_F 1.0
constexpr real1_f PI_R1 = (real1_f)M_PI;
constexpr real1_f SQRT2_R1 = (real1_f)M_SQRT2;
constexpr real1_f SQRT1_2_R1 = (real1_f)M_SQRT1_2;
#define REAL1_DEFAULT_ARG -999.0
// Half the probability in any single permutation of 192 maximally superposed qubits
#define REAL1_EPSILON 7.965459555662261e-59
#endif

#if ENABLE_CUDA
#if FPPOW < 5
#include <cuda_fp16.h>
#define qCudaReal1 __half
#define qCudaReal2 __half2
#define qCudaReal4 __half2*
#define qCudaCmplx __half2
#define qCudaCmplx2 __half2*
#define qCudaReal1_f float
#define make_qCudaCmplx make_half2
#define ZERO_R1_CUDA ((qCudaReal1)0.0f)
#define REAL1_EPSILON_CUDA 0.000000477f
#define PI_R1_CUDA M_PI
#elif FPPOW < 6
#define qCudaReal1 float
#define qCudaReal2 float2
#define qCudaReal4 float4
#define qCudaCmplx float2
#define qCudaCmplx2 float4
#define qCudaReal1_f float
#define make_qCudaCmplx make_float2
#define make_qCudaCmplx2 make_float4
#define ZERO_R1_CUDA 0.0f
#define REAL1_EPSILON_CUDA REAL1_EPSILON
#define PI_R1_CUDA PI_R1
#else
#define qCudaReal1 double
#define qCudaReal2 double2
#define qCudaReal4 double4
#define qCudaCmplx double2
#define qCudaCmplx2 double4
#define qCudaReal1_f double
#define make_qCudaCmplx make_double2
#define make_qCudaCmplx2 make_double4
#define ZERO_R1_CUDA 0.0
#define REAL1_EPSILON_CUDA REAL1_EPSILON
#define PI_R1_CUDA PI_R1
#endif
#endif

constexpr complex ONE_CMPLX = { ONE_R1, ZERO_R1 };
constexpr complex ZERO_CMPLX = { ZERO_R1, ZERO_R1 };
constexpr complex I_CMPLX = { ZERO_R1, ONE_R1 };
constexpr complex CMPLX_DEFAULT_ARG = { REAL1_DEFAULT_ARG, REAL1_DEFAULT_ARG };

constexpr real1 FP_NORM_EPSILON = (real1)(std::numeric_limits<real1>::epsilon() / 2);
constexpr real1_f TRYDECOMPOSE_EPSILON = (real1_f)(16 * FP_NORM_EPSILON);
constexpr real1_f FP_NORM_EPSILON_F = std::numeric_limits<real1_f>::epsilon() / 2;

} // namespace Qrack

