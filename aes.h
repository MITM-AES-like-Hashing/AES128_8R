#ifndef AES_H__
#define AES_H__
#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include <set>
#include <array>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <immintrin.h>

using namespace std;

typedef unsigned char   word8;
typedef unsigned long   word32;

#if defined(_MSC_VER)
#define ALIGNED_(x) __declspec(align(x))
#else
#if defined(__GNUC__)
#define ALIGNED_(x) __attribute__((aligned(x)))
typedef unsigned long long _ULonglong;
#endif
#endif

#define ALIGNED_TYPE_(t,x) t ALIGNED_(x)


word8 Logtable[256] = {
  0,   0,  25,   1,  50,   2,  26, 198,  75, 199,  27, 104,  51, 238, 223,   3, 
100,   4, 224,  14,  52, 141, 129, 239,  76, 113,   8, 200, 248, 105,  28, 193, 
125, 194,  29, 181, 249, 185,  39, 106,  77, 228, 166, 114, 154, 201,   9, 120, 
101,  47, 138,   5,  33,  15, 225,  36,  18, 240, 130,  69,  53, 147, 218, 142, 
150, 143, 219, 189,  54, 208, 206, 148,  19,  92, 210, 241,  64,  70, 131,  56, 
102, 221, 253,  48, 191,   6, 139,  98, 179,  37, 226, 152,  34, 136, 145,  16, 
126, 110,  72, 195, 163, 182,  30,  66,  58, 107,  40,  84, 250, 133,  61, 186, 
 43, 121,  10,  21, 155, 159,  94, 202,  78, 212, 172, 229, 243, 115, 167,  87, 
175,  88, 168,  80, 244, 234, 214, 116,  79, 174, 233, 213, 231, 230, 173, 232, 
 44, 215, 117, 122, 235,  22,  11, 245,  89, 203,  95, 176, 156, 169,  81, 160, 
127,  12, 246, 111,  23, 196,  73, 236, 216,  67,  31,  45, 164, 118, 123, 183, 
204, 187,  62,  90, 251,  96, 177, 134,  59,  82, 161, 108, 170,  85,  41, 157, 
151, 178, 135, 144,  97, 190, 220, 252, 188, 149, 207, 205,  55,  63,  91, 209, 
 83,  57, 132,  60,  65, 162, 109,  71,  20,  42, 158,  93,  86, 242, 211, 171, 
 68,  17, 146, 217,  35,  32,  46, 137, 180, 124, 184,  38, 119, 153, 227, 165, 
103,  74, 237, 222, 197,  49, 254,  24,  13,  99, 140, 128, 192, 247, 112,   7, 
};

word8 Alogtable[256] = {
  1,   3,   5,  15,  17,  51,  85, 255,  26,  46, 114, 150, 161, 248,  19,  53, 
 95, 225,  56,  72, 216, 115, 149, 164, 247,   2,   6,  10,  30,  34, 102, 170, 
229,  52,  92, 228,  55,  89, 235,  38, 106, 190, 217, 112, 144, 171, 230,  49, 
 83, 245,   4,  12,  20,  60,  68, 204,  79, 209, 104, 184, 211, 110, 178, 205, 
 76, 212, 103, 169, 224,  59,  77, 215,  98, 166, 241,   8,  24,  40, 120, 136, 
131, 158, 185, 208, 107, 189, 220, 127, 129, 152, 179, 206,  73, 219, 118, 154, 
181, 196,  87, 249,  16,  48,  80, 240,  11,  29,  39, 105, 187, 214,  97, 163, 
254,  25,  43, 125, 135, 146, 173, 236,  47, 113, 147, 174, 233,  32,  96, 160, 
251,  22,  58,  78, 210, 109, 183, 194,  93, 231,  50,  86, 250,  21,  63,  65, 
195,  94, 226,  61,  71, 201,  64, 192,  91, 237,  44, 116, 156, 191, 218, 117, 
159, 186, 213, 100, 172, 239,  42, 126, 130, 157, 188, 223, 122, 142, 137, 128, 
155, 182, 193,  88, 232,  35, 101, 175, 234,  37, 111, 177, 200,  67, 197,  84, 
252,  31,  33,  99, 165, 244,   7,   9,  27,  45, 119, 153, 176, 203,  70, 202, 
 69, 207,  74, 222, 121, 139, 134, 145, 168, 227,  62,  66, 198,  81, 243,  14, 
 18,  54,  90, 238,  41, 123, 141, 140, 143, 138, 133, 148, 167, 242,  13,  23, 
 57,  75, 221, 124, 132, 151, 162, 253,  28,  36, 108, 180, 199,  82, 246,   1, 
};

word8 S[256] = {
 99, 124, 119, 123, 242, 107, 111, 197,  48,   1, 103,  43, 254, 215, 171, 118, 
202, 130, 201, 125, 250,  89,  71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 
183, 253, 147,  38,  54,  63, 247, 204,  52, 165, 229, 241, 113, 216,  49,  21, 
  4, 199,  35, 195,  24, 150,   5, 154,   7,  18, 128, 226, 235,  39, 178, 117, 
  9, 131,  44,  26,  27, 110,  90, 160,  82,  59, 214, 179,  41, 227,  47, 132, 
 83, 209,   0, 237,  32, 252, 177,  91, 106, 203, 190,  57,  74,  76,  88, 207, 
208, 239, 170, 251,  67,  77,  51, 133,  69, 249,   2, 127,  80,  60, 159, 168, 
 81, 163,  64, 143, 146, 157,  56, 245, 188, 182, 218,  33,  16, 255, 243, 210, 
205,  12,  19, 236,  95, 151,  68,  23, 196, 167, 126,  61, 100,  93,  25, 115, 
 96, 129,  79, 220,  34,  42, 144, 136,  70, 238, 184,  20, 222,  94,  11, 219, 
224,  50,  58,  10,  73,   6,  36,  92, 194, 211, 172,  98, 145, 149, 228, 121, 
231, 200,  55, 109, 141, 213,  78, 169, 108,  86, 244, 234, 101, 122, 174,   8, 
186, 120,  37,  46,  28, 166, 180, 198, 232, 221, 116,  31,  75, 189, 139, 138, 
112,  62, 181, 102,  72,   3, 246,  14,  97,  53,  87, 185, 134, 193,  29, 158, 
225, 248, 152,  17, 105, 217, 142, 148, 155,  30, 135, 233, 206,  85,  40, 223, 
140, 161, 137,  13, 191, 230,  66, 104,  65, 153,  45,  15, 176,  84, 187,  22, 
};

word8 iS[256] = {
 82,   9, 106, 213,  48,  54, 165,  56, 191,  64, 163, 158, 129, 243, 215, 251, 
124, 227,  57, 130, 155,  47, 255, 135,  52, 142,  67,  68, 196, 222, 233, 203, 
 84, 123, 148,  50, 166, 194,  35,  61, 238,  76, 149,  11,  66, 250, 195,  78, 
  8,  46, 161, 102,  40, 217,  36, 178, 118,  91, 162,  73, 109, 139, 209,  37, 
114, 248, 246, 100, 134, 104, 152,  22, 212, 164,  92, 204,  93, 101, 182, 146, 
108, 112,  72,  80, 253, 237, 185, 218,  94,  21,  70,  87, 167, 141, 157, 132, 
144, 216, 171,   0, 140, 188, 211,  10, 247, 228,  88,   5, 184, 179,  69,   6, 
208,  44,  30, 143, 202,  63,  15,   2, 193, 175, 189,   3,   1,  19, 138, 107, 
 58, 145,  17,  65,  79, 103, 220, 234, 151, 242, 207, 206, 240, 180, 230, 115, 
150, 172, 116,  34, 231, 173,  53, 133, 226, 249,  55, 232,  28, 117, 223, 110, 
 71, 241,  26, 113,  29,  41, 197, 137, 111, 183,  98,  14, 170,  24, 190,  27, 
252,  86,  62,  75, 198, 210, 121,  32, 154, 219, 192, 254, 120, 205,  90, 244, 
 31, 221, 168,  51, 136,   7, 199,  49, 177,  18,  16,  89,  39, 128, 236,  95, 
 96,  81, 127, 169,  25, 181,  74,  13,  45, 229, 122, 159, 147, 201, 156, 239, 
160, 224,  59,  77, 174,  42, 245, 176, 200, 235, 187,  60, 131,  83, 153,  97, 
 23,  43,   4, 126, 186, 119, 214,  38, 225, 105,  20,  99,  85,  33,  12, 125, 
};

word8 iG[4][4] = {
0x0e, 0x09, 0x0d, 0x0b, 
0x0b, 0x0e, 0x09, 0x0d, 
0x0d, 0x0b, 0x0e, 0x09, 
0x09, 0x0d, 0x0b, 0x0e, 
};

word32 rcon[30] = { 
  0x01,0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a, 0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91, };

word8 mul(word8 a, word8 b) {
   /* multiply two elements of GF(2^m)
    * needed for MixColumn and InvMixColumn
    */
	if (a && b) return Alogtable[(Logtable[a] + Logtable[b])%255];
	else return 0;
}

void mix_column(array<unsigned char, 4> &out, array<unsigned char, 4> &in)
{
    array<unsigned char, 4> b;
    unsigned char c;
    unsigned char h;
    /* The array 'a' is simply a copy of the input array 'r'
     * The array 'b' is each element of the array 'a' multiplied by 2
     * in Rijndael's Galois field
     * a[n] ^ b[n] is element n multiplied by 3 in Rijndael's Galois field */ 
    for (c = 0; c < 4; c++) {
        /* h is 0xff if the high bit of in[c] is set, 0 otherwise */
        h = (unsigned char)((signed char)in[c] >> 7); /* arithmetic right shift, thus shifting in either zeros or ones */
        b[c] = in[c] << 1; /* implicitly removes high bit because b[c] is an 8-bit char, so we xor by 0x1b and not 0x11b in the next line */
        b[c] ^= 0x1B & h; /* Rijndael's Galois field */
    }
    out[0] = b[0] ^ in[3] ^ in[2] ^ b[1] ^ in[1]; /* 2 * a0 + a3 + a2 + 3 * a1 */
    out[1] = b[1] ^ in[0] ^ in[3] ^ b[2] ^ in[2]; /* 2 * a1 + a0 + a3 + 3 * a2 */
    out[2] = b[2] ^ in[1] ^ in[0] ^ b[3] ^ in[3]; /* 2 * a2 + a1 + a0 + 3 * a3 */
    out[3] = b[3] ^ in[2] ^ in[1] ^ b[0] ^ in[0]; /* 2 * a3 + a2 + a1 + 3 * a0 */
}

inline __m128i mix_column(__m128i &in)
{
    __m128i out;
    out = _mm_aesdeclast_si128(in, _mm_setzero_si128());
    out = _mm_aesenc_si128(out, _mm_setzero_si128());
    return out;
}

inline __m128i inv_mix_column(__m128i &in)
{
    return _mm_aesimc_si128(in);
}

void inv_mix_column(array<unsigned char, 4> &out, array<unsigned char, 4> &in) {
    out[0] = mul(0xe, in[0]) ^ mul(0xb, in[1]) ^ mul(0xd, in[2]) ^ mul(0x9, in[3]);
    out[1] = mul(0xe, in[1]) ^ mul(0xb, in[2]) ^ mul(0xd, in[3]) ^ mul(0x9, in[0]);
    out[2] = mul(0xe, in[2]) ^ mul(0xb, in[3]) ^ mul(0xd, in[0]) ^ mul(0x9, in[1]);
    out[3] = mul(0xe, in[3]) ^ mul(0xb, in[0]) ^ mul(0xd, in[1]) ^ mul(0x9, in[2]);
}

unsigned char inv_mix_column_cell0(array<unsigned char, 4> &in) {
    return mul(0xe, in[0]) ^ mul(0xb, in[1]) ^ mul(0xd, in[2]) ^ mul(0x9, in[3]);
}
unsigned char inv_mix_column_cell1(array<unsigned char, 4> &in) {
    return mul(0xe, in[1]) ^ mul(0xb, in[2]) ^ mul(0xd, in[3]) ^ mul(0x9, in[0]);
}
unsigned char inv_mix_column_cell2(array<unsigned char, 4> &in) {
    return mul(0xe, in[2]) ^ mul(0xb, in[3]) ^ mul(0xd, in[0]) ^ mul(0x9, in[1]);
}
unsigned char inv_mix_column_cell3(array<unsigned char, 4> &in) {
    return mul(0xe, in[3]) ^ mul(0xb, in[0]) ^ mul(0xd, in[1]) ^ mul(0x9, in[2]);
}

#define reverse(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xa,xb,xc,xd,xe,xf) xf,xe,xd,xc,xb,xa,x9,x8,x7,x6,x5,x4,x3,x2,x1,x0

const __m128i RC[] =
{
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x01),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x02),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x04),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x08),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x10),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x20),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x40),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x80),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x1b),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x36),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x6c),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xd8),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xab),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x4d),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x9a),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x2f),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x5e),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xbc),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x63),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xc6),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x97),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x35),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x6a),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xd4),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xb3),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x7d),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xfa),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xef),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0xc5),
  _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0x91),
};

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

const __m128i KEY_SHUFFLE_OFFSET = _mm_set_epi8(State_to_Parameter_Order(
                                                  13, 0xFF, 0xFF, 0xFF,
                                                0xFF,   14, 0xFF, 0xFF,
                                                0xFF, 0xFF,   15, 0xFF,
                                                0xFF, 0xFF, 0xFF,   12
                                                )); //  12, 0xFF, 0xFF, 0xFF, 0xFF,   15, 0xFF, 0xFF, 0xFF, 0xFF,   14, 0xFF, 0xFF, 0xFF, 0xFF,   13
const __m128i zeroKey = _mm_set_epi8(State_to_Parameter_Order(
                                     _ , 0x63, 0x63, 0x63, 
                                     _ , 0x63, 0x63, 0x63,
                                     _ , 0x63, 0x63, 0x63, 
                                     _ , 0x63, 0x63, 0x63
                                     )); // 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x63, 0x63, 0x63, 0x63

__m128i AES128_key_exp(__m128i &in, int r)
{
    // 0, 4, 8,12
    // 1, 5, 9,13
    // 2, 6,10,14
    // 3, 7,11,15
    //
    //13, 0, 4, 8,
    //14, 1, 5, 9,
    //15, 2, 6,10,
    //12, 3, 7,11,
    //
    //13,FF,FF,FF,
    //FF,14,FF,FF,
    //FF,FF,15,FF,
    //FF,FF,FF,12,
    __m128i out;

    out = _mm_shuffle_epi8(in, KEY_SHUFFLE_OFFSET);
    /* out = 
     * in[13],      0,      0,     0,
     *      0, in[14],      0,     0,
     *      0,      0, in[15],     0,
     *      0,      0,      0, in[12]
    */
    out = _mm_aesenclast_si128(out, zeroKey);
    /* out = 
     * S[in[13],      0,      0,     0, 
     * S[in[14],      0,      0,     0, 
     * S[in[15],      0,      0,     0, 
     * S[in[12],      0,      0,     0, 
    */
    out = _mm_xor_si128(out, in);
    out = _mm_xor_si128(out, RC[r]);
    out = _mm_xor_si128(out, _mm_slli_si128(out, 4)); // col0 = col0, col1 = col1 ^ col0, col2 = col2 ^ col1;        col3 = col3 ^ col2
    out = _mm_xor_si128(out, _mm_slli_si128(out, 4)); // col0 = col0, col1 = col1,        col2 = col2 ^ col0;        col3 = col3 ^ col1
    out = _mm_xor_si128(out, _mm_slli_si128(out, 4)); // col0 = col0, col1 = col1 ^ col0, col2 = col2 ^ col0 ^ col1; col3 = col3 ^ col1 ^ col2 ^ col0
    return out;
}

__m128i invAES128_key_exp(__m128i &in, int r)
{
    __m128i out;
    __m128i tmp;
    out = _mm_xor_si128(in, _mm_slli_si128(in, 4)); // col0 = col0, col1 = col1 ^ col0, col2 = col2 ^ col1;        col3 = col3 ^ col2
    out = _mm_xor_si128(out, RC[r]);
    tmp = _mm_shuffle_epi8(out, KEY_SHUFFLE_OFFSET);
    tmp = _mm_aesenclast_si128(tmp, zeroKey);
    out = _mm_xor_si128(out, tmp);
    return out;
}

// Another way of AES KeySchedule, using _mm_aeskeygenassist_si128
// which is slower
#define AES_128_key_exp(k, rconi) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rconi))
static __m128i aes_128_key_expansion(__m128i key, __m128i keygened){
    keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
    return _mm_xor_si128(key, keygened);
}
//public API
static void aes128_load_key_and_exp(uint8_t *key, __m128i *rnkey){
    rnkey[0] = _mm_loadu_si128((const __m128i*) key);
    rnkey[1]  = AES_128_key_exp(rnkey[0], 0x01);
    rnkey[2]  = AES_128_key_exp(rnkey[1], 0x02);
    rnkey[3]  = AES_128_key_exp(rnkey[2], 0x04);
    rnkey[4]  = AES_128_key_exp(rnkey[3], 0x08);
    rnkey[5]  = AES_128_key_exp(rnkey[4], 0x10);
    rnkey[6]  = AES_128_key_exp(rnkey[5], 0x20);
    rnkey[7]  = AES_128_key_exp(rnkey[6], 0x40);
    rnkey[8]  = AES_128_key_exp(rnkey[7], 0x80);
    rnkey[9]  = AES_128_key_exp(rnkey[8], 0x1B);
    rnkey[10] = AES_128_key_exp(rnkey[9], 0x36);
}

#endif