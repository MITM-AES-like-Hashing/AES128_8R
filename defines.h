#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdint.h>
#include <immintrin.h>
#include <array>
#include <tuple>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <csetjmp>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "omp.h"

using namespace std;

// Cipher Parameters
#define ROW_N (4)
#define COL_N (4)
#define STATE_N (ROW_N * COL_N)

#define SB_n (8U)
#define SB_a (1U << (SB_n))

// Attack configurations
#define k4_CONST_SAMPLE_N (1U << (3 * 8))
#define MC3_CONST_SAMPLE_N (1U << (12 * 8))
#define MC2_CONST_SAMPLE_N (1U << (3 * 8))
#define AK5_CONST_SAMPLE_N (1U << (2 * 8))
#define AK4_CONST_COL_SAMPLE_N (1U << (3 * 8))

#ifndef DoFF
#define DoFF 8
#endif

#if (DoFF == 16)
#define INIT_k3_col0_N (SB_a)
#define INIT_k3_col1_N (SB_a)
#define INIT_k3_col2_N (SB_a)
#define INIT_k3_col3_N (SB_a)
#define INIT_SB4_cell00_N (SB_a)
#define INIT_SB4_cell05_N (SB_a)
#define INIT_SB4_cell10_N (SB_a)
#define INIT_SB4_cell15_N (SB_a)
#elif (DoFF == 8)
#define INIT_k3_col0_N (SB_a)
#define INIT_k3_col1_N (SB_a)
#define INIT_k3_col2_N (SB_a)
#define INIT_k3_col3_N (SB_a)
#define INIT_SB4_cell00_N (SB_a / (1U << 2))
#define INIT_SB4_cell05_N (SB_a / (1U << 2))
#define INIT_SB4_cell10_N (SB_a / (1U << 2))
#define INIT_SB4_cell15_N (SB_a / (1U << 2))
#else // (DoFF == 4)
#define INIT_k3_col0_N (SB_a)
#define INIT_k3_col1_N (SB_a)
#define INIT_k3_col2_N (SB_a)
#define INIT_k3_col3_N (SB_a)
#define INIT_SB4_cell00_N (SB_a / (1U << 3))
#define INIT_SB4_cell05_N (SB_a / (1U << 3))
#define INIT_SB4_cell10_N (SB_a / (1U << 3))
#define INIT_SB4_cell15_N (SB_a / (1U << 3))
#endif

#define INIT_k3_col0_SB4_cell00_N ((INIT_k3_col0_N) * (INIT_SB4_cell00_N))
#define INIT_k3_col1_SB4_cell05_N ((INIT_k3_col1_N) * (INIT_SB4_cell05_N))
#define INIT_k3_col2_SB4_cell10_N ((INIT_k3_col2_N) * (INIT_SB4_cell10_N))
#define INIT_k3_col3_SB4_cell15_N ((INIT_k3_col3_N) * (INIT_SB4_cell15_N))


const __m128i one_128 = _mm_set_epi64x(0xffffffffffffffffULL, 0xffffffffffffffffULL);

#ifndef State_to_Parameter_Order
// 00, 04, 08, 12,
// 01, 05, 09, 13,
// 02, 06, 10, 14,
// 03, 07, 11, 15,
#define State_to_Parameter_Order( \
    x00, x04, x08, x12, \
    x01, x05, x09, x13, \
    x02, x06, x10, x14, \
    x03, x07, x11, x15) \
    x15, x14, x13, x12, x11, x10, x09, x08, x07, x06, x05, x04, x03, x02, x01, x00
#endif

#ifndef _
#define _ (0x00)
#endif

#ifndef PARTIAL
#define PARTIAL 32
#endif

#if (PARTIAL == 64)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ , 0xFF,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0xFF,   _ , 0xFF,   _ ,
                                                        0xFF,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL == 56)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0xFF,   _ , 0xFF,   _ ,
                                                        0xFF,   _ ,   _ , 0xFF
                                                        ));
#elif (PARTIAL == 48)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0xFF,   _ , 0xFF,   _ ,
                                                        0xFF,   _ ,   _ ,   _ 
                                                        ));
#elif (PARTIAL == 40)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF, 0xFF,   _ ,   _ ,
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF,   _ ,   _ ,   _ 
                                                        ));
#else //(PARTIAL == 32)
const __m128i PARTIAL_MATCH_TARGET_MASK = _mm_set_epi8(State_to_Parameter_Order(
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF,   _ ,   _ ,   _ ,
                                                        0xFF,   _ ,   _ ,   _ 
                                                        ));
#endif

#define PRINT_COUNT 0
#define DEBUG 0

#ifndef omp_nb_threads
#define omp_nb_threads (4)
#endif

typedef uint8_t cell_t;
typedef array<cell_t, ROW_N> column_t;
typedef __m128i state_t;

typedef struct
{
    column_t k3_col;
    cell_t MC3_W;
    cell_t SB4_B;
} colcell_t;

template<typename T>
inline T operator ^ (const T& lhs, const T& rhs)
{
  T tmp;
  for (size_t i = 0; i < lhs.size(); i++)
  {
    tmp[i] = lhs[i] ^ rhs[i];
  }
  return tmp;
}

#define EQU(w, u)  (_mm_testz_si128(one_128, _mm_xor_si128(w, u)) == 1)
#define GE(w, u)  (EQU(one_128, _mm_cmpgt_epi32(w, u)))

struct cmpState_t {
    bool operator()(const state_t& a, const state_t& b) const {
        return GE(a, b);
    }
};

struct cmpColcell_t {
    bool operator()(const colcell_t& a, const colcell_t& b) const {
        return (memcmp(&a, &b, sizeof(colcell_t)) > 0);
    }
};

inline void printstate(state_t &state)
{
    cell_t * state_pt = (cell_t *)&state;

    cout << hex << setfill('0');
    for (size_t ri = 0; ri < ROW_N; ri++)
    {
        for (size_t ci = 0; ci < COL_N; ci++)
        {

            cout << setw(sizeof(cell_t) * 2) << state_pt[ci * ROW_N + ri] + '\0' << "  ";
        }
        cout << endl;
    }
    cout << endl;
    cout << dec << setfill(' ');
}