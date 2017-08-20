/*
 * Tiny self-contained version of the PCG Random Number Generation for C++ put
 * together from pieces of the much larger C/C++ codebase with vectorization
 * using Enoki.
 *
 * Wenzel Jakob, February 2017
 *
 * The PCG random number generator was developed by Melissa O'Neill
 * <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */

#pragma once

#include <mitsuba/core/simd.h>
#include <mitsuba/core/logger.h>
#include <enoki/random.h>

NAMESPACE_BEGIN(mitsuba)

using PCG32  = enoki::PCG32<uint32_t>;
using PCG32P = enoki::PCG32<UInt32P>;

/**
 * \brief Generate fast and reasonably good pseudorandom numbers using the
 * Tiny Encryption Algorithm (TEA) by David Wheeler and Roger Needham.
 *
 * For details, refer to "GPU Random Numbers via the Tiny Encryption Algorithm"
 * by Fahad Zafar, Marc Olano, and Aaron Curtis.
 *
 * \param v0
 *     First input value to be encrypted (could be the sample index)
 * \param v1
 *     Second input value to be encrypted (e.g. the requested random number dimension)
 * \param rounds
 *     How many rounds should be executed? The default for random number
 *     generation is 4.
 * \return
 *     A uniformly distributed 32-bit integer
 */

template <typename UInt32>
UInt32 sample_tea_32(UInt32 v0, UInt32 v1, int rounds = 4) {
    UInt32 sum = 0;

    ENOKI_NOUNROLL for (int i = 0; i < rounds; ++i) {
        sum += 0x9e3779b9;
        v0 += (sli<4>(v1) + 0xa341316c) ^ (v1 + sum) ^ (sri<5>(v1) + 0xc8013ea4);
        v1 += (sli<4>(v0) + 0xad90777d) ^ (v0 + sum) ^ (sri<5>(v0) + 0x7e95761e);
    }

    return v1;
}

extern template MTS_EXPORT_CORE uint32_t sample_tea_32(uint32_t, uint32_t, int);
extern template MTS_EXPORT_CORE UInt32P  sample_tea_32(UInt32P,  UInt32P,  int);

/**
 * \brief Generate fast and reasonably good pseudorandom numbers using the
 * Tiny Encryption Algorithm (TEA) by David Wheeler and Roger Needham.
 *
 * For details, refer to "GPU Random Numbers via the Tiny Encryption Algorithm"
 * by Fahad Zafar, Marc Olano, and Aaron Curtis.
 *
 * \param v0
 *     First input value to be encrypted (could be the sample index)
 * \param v1
 *     Second input value to be encrypted (e.g. the requested random number dimension)
 * \param rounds
 *     How many rounds should be executed? The default for random number
 *     generation is 4.
 * \return
 *     A uniformly distributed 64-bit integer
 */

template <typename UInt32>
uint64_array_t<UInt32> sample_tea_64(UInt32 v0, UInt32 v1, int rounds = 4) {
    UInt32 sum = 0;

    ENOKI_NOUNROLL for (int i = 0; i < rounds; ++i) {
        sum += 0x9e3779b9;
        v0 += (sli<4>(v1) + 0xa341316c) ^ (v1 + sum) ^ (sri<5>(v1) + 0xc8013ea4);
        v1 += (sli<4>(v0) + 0xad90777d) ^ (v0 + sum) ^ (sri<5>(v0) + 0x7e95761e);
    }

    return uint64_array_t<UInt32>(v0) + sli<32>(uint64_array_t<UInt32>(v1));
}

extern template MTS_EXPORT_CORE uint64_array_t<uint32_t> sample_tea_64(uint32_t, uint32_t, int);
extern template MTS_EXPORT_CORE uint64_array_t<UInt32P>  sample_tea_64(UInt32P,  UInt32P,  int);

/**
 * \brief Generate fast and reasonably good pseudorandom numbers using the
 * Tiny Encryption Algorithm (TEA) by David Wheeler and Roger Needham.
 *
 * This function uses \ref sample_tea_ to return single precision floating point
 * numbers on the interval <tt>[0, 1)</tt>
 *
 * \param v0
 *     First input value to be encrypted (could be the sample index)
 * \param v1
 *     Second input value to be encrypted (e.g. the requested random number dimension)
 * \param rounds
 *     How many rounds should be executed? The default for random number
 *     generation is 4.
 * \return
 *     A uniformly distributed floating point number on the interval <tt>[0, 1)</tt>
 */
template <typename UInt32>
float32_array_t<UInt32> sample_tea_float32(UInt32 v0, UInt32 v1, int rounds = 4) {
    return reinterpret_array<float32_array_t<UInt32>>(
        sri<9>(sample_tea_32(v0, v1, rounds)) | 0x3f800000u) - 1.f;
}

extern template MTS_EXPORT_CORE float32_array_t<uint32_t> sample_tea_float32(uint32_t, uint32_t, int);
extern template MTS_EXPORT_CORE float32_array_t<UInt32P>  sample_tea_float32(UInt32P,  UInt32P,  int);

/**
 * \brief Generate fast and reasonably good pseudorandom numbers using the
 * Tiny Encryption Algorithm (TEA) by David Wheeler and Roger Needham.
 *
 * This function uses \ref sample_tea_ to return double precision floating point
 * numbers on the interval <tt>[0, 1)</tt>
 *
 * \param v0
 *     First input value to be encrypted (could be the sample index)
 * \param v1
 *     Second input value to be encrypted (e.g. the requested random number dimension)
 * \param rounds
 *     How many rounds should be executed? The default for random number
 *     generation is 4.
 * \return
 *     A uniformly distributed floating point number on the interval <tt>[0, 1)</tt>
 */

template <typename UInt32>
float64_array_t<UInt32> sample_tea_float64(UInt32 v0, UInt32 v1, int rounds = 4) {
    return reinterpret_array<float64_array_t<UInt32>>(
        sri<12>(sample_tea_64(v0, v1, rounds)) | 0x3ff0000000000000ull) - 1.0;
}

extern template MTS_EXPORT_CORE float64_array_t<uint32_t> sample_tea_float64(uint32_t, uint32_t, int);
extern template MTS_EXPORT_CORE float64_array_t<UInt32P>  sample_tea_float64(UInt32P,  UInt32P,  int);

/// Alias to \ref sample_tea_float32 or \ref sample_tea_float64 based on compilation flags
template <typename UInt32>
auto sample_tea_float(UInt32 v0, UInt32 v1, int rounds = 4) {
    #if defined(SINGLE_PRECISION)
        return sample_tea_float32(v0, v1, rounds);
    #else
        return sample_tea_float64(v0, v1, rounds);
    #endif
}

NAMESPACE_END(mitsuba)
