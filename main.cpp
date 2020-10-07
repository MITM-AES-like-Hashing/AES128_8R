#include "rndfill.h"
#include "aes.h"
#include "defines.h"

using namespace std;

array<vector<uint64_t>, k4_CONST_SAMPLE_N> k4_C_k3SB4_B;

/**************************************************************************************
 * Precompute lookup table k4_C_k3SB4_B,
 * the index is the values of the 3 Gray bytes in k4, 
 * the value is the Blue cells (neutral bytes for forward) in state #SB4 and
 * the bytes to generate Blue bytes in k3.
 * To save memory (only by a constant factor), we do not store k4 directly.
 **************************************************************************************/
void precomputeB(
    array<uint8_t, 12> &MC3_C,
    array<uint8_t,  3> &MC2_C)
{
    cout << "Start precompute the Blue cells (neutral bytes for forward) in state k4 and #SB4." << endl;
    uint64_t complexity_tmp = 0ULL;
    uint64_t complexity = 0ULL;

    // Precompute BLUE
    column_t tmp = {0};
    column_t MC3_00_col0_MCout;
    column_t MC3_07_col1_MCout;
    column_t MC3_10_col2_MCout;
    column_t MC3_13_col3_MCout;

    array<colcell_t, INIT_k3_col0_SB4_cell00_N > k3_col0_SB4_cell00;
    array<colcell_t, INIT_k3_col1_SB4_cell05_N > k3_col1_SB4_cell05;
    array<colcell_t, INIT_k3_col2_SB4_cell10_N > k3_col2_SB4_cell10;
    array<colcell_t, INIT_k3_col3_SB4_cell15_N > k3_col3_SB4_cell15;

    multimap<colcell_t, colcell_t, cmpColcell_t> k3_col1_0;
    multimap<colcell_t, colcell_t, cmpColcell_t> k3_col2_1;
    multimap<colcell_t, colcell_t, cmpColcell_t> k3_col2_3;

    multimap<colcell_t, pair<colcell_t, colcell_t >, cmpColcell_t> k3_col1_23;

    typedef std::multimap<colcell_t, colcell_t, cmpColcell_t>::iterator MAP2_3_Iterator;
    typedef std::multimap<colcell_t, pair<colcell_t, colcell_t >, cmpColcell_t>::iterator MAP1_23_Iterator;
    colcell_t colcell_tmp;

    // invMixColumns((#SB4[ 0],        0,        0,        0) ⊕ (k3[ 0], k3[ 1], k3[ 2], k3[ 3])) = (        *, MC3_C[ 0], MC3_C[ 1], MC3_C[ 2])
    // invMixColumns((       0, #SB4[ 5],        0,        0) ⊕ (k3[ 4], k3[ 5], k3[ 6], k3[ 7])) = (MC3_C[ 3], MC3_C[ 4], MC3_C[ 5],         *)
    // invMixColumns((       0,        0, #SB4[10],        0) ⊕ (k3[ 8], k3[ 9], k3[10], k3[11])) = (MC3_C[ 6], MC3_C[ 7],         *, MC3_C[ 8])
    // invMixColumns((       0,        0,        0, #SB4[15]) ⊕ (k3[12], k3[13], k3[14], k3[15])) = (MC3_C[ 9],         *, MC3_C[10], MC3_C[11])

    // (#SB4[ 0],        0,        0,        0) ⊕ (k3[ 0], k3[ 1], k3[ 2], k3[ 3]) = MixColumns(        *, MC3_C[ 0], MC3_C[ 1], MC3_C[ 2])
    // (       0, #SB4[ 5],        0,        0) ⊕ (k3[ 4], k3[ 5], k3[ 6], k3[ 7]) = MixColumns(MC3_C[ 3], MC3_C[ 4], MC3_C[ 5],         *)
    // (       0,        0, #SB4[10],        0) ⊕ (k3[ 8], k3[ 9], k3[10], k3[11]) = MixColumns(MC3_C[ 6], MC3_C[ 7],         *, MC3_C[ 8])
    // (       0,        0,        0, #SB4[15]) ⊕ (k3[12], k3[13], k3[14], k3[15]) = MixColumns(MC3_C[ 9],         *, MC3_C[10], MC3_C[11])

    // (k3[ 0], k3[ 1], k3[ 2], k3[ 3]) = MixColumns(        *, MC3_C[ 0], MC3_C[ 1], MC3_C[ 2]) ⊕ (#SB4[ 0],        0,        0,        0)
    // (k3[ 4], k3[ 5], k3[ 6], k3[ 7]) = MixColumns(MC3_C[ 3], MC3_C[ 4], MC3_C[ 5],         *) ⊕ (       0, #SB4[ 5],        0,        0)
    // (k3[ 8], k3[ 9], k3[10], k3[11]) = MixColumns(MC3_C[ 6], MC3_C[ 7],         *, MC3_C[ 8]) ⊕ (       0,        0, #SB4[10],        0)
    // (k3[12], k3[13], k3[14], k3[15]) = MixColumns(MC3_C[ 9],         *, MC3_C[10], MC3_C[11]) ⊕ (       0,        0,        0, #SB4[15])
    complexity_tmp = 0ULL;
    for (size_t i = 0; i < INIT_k3_col0_N; i++)
    {
        // enumerate on cell 0 in state #MC3
        // (        *, MC3_C[ 0], MC3_C[ 1], MC3_C[ 2])
        tmp[0] = (uint8_t) i & 0xffU;
        tmp[1] = MC3_C[0];
        tmp[2] = MC3_C[1];
        tmp[3] = MC3_C[2];
        mix_column(MC3_00_col0_MCout, tmp);

        colcell_tmp.MC3_W = tmp[0];
        colcell_tmp.k3_col = MC3_00_col0_MCout;
        colcell_tmp.SB4_B = 0U;
        k3_col0_SB4_cell00[i * INIT_SB4_cell00_N] = colcell_tmp;
        for (size_t j = 1; j < INIT_SB4_cell00_N; j++)
        {
            // enumerate on cell 0 in state #SB4
            // obtain candidates for column 0 of k3 and cell 0 of #SB4
            MC3_00_col0_MCout[0] ^= ((cell_t)(j - 1) ^ (cell_t)j);
            colcell_tmp.k3_col = MC3_00_col0_MCout;
            colcell_tmp.SB4_B = (uint8_t) j & 0xffU;
            k3_col0_SB4_cell00[i * INIT_SB4_cell00_N + j] = colcell_tmp;
            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;

    complexity_tmp = 0ULL;
    for (size_t i = 0; i < INIT_k3_col1_N; i++)
    {
        // enumerate on cell 7 in state #MC3
        // (MC3_C[ 3], MC3_C[ 4], MC3_C[ 5],         *)
        tmp[0] = MC3_C[3];
        tmp[1] = MC3_C[4];
        tmp[2] = MC3_C[5];
        tmp[3] = (uint8_t) i & 0xffU;
        mix_column(MC3_07_col1_MCout, tmp);

        colcell_tmp.MC3_W = tmp[3];
        colcell_tmp.k3_col = MC3_07_col1_MCout;
        colcell_tmp.SB4_B = 0U;
        k3_col1_SB4_cell05[i * INIT_SB4_cell05_N] = colcell_tmp;
        for (size_t j = 1; j < INIT_SB4_cell05_N; j++)
        {
            // enumerate on cell 5 in state #SB4
            // obtain candidates for column 1 of k3 and cell 5 of #SB4
            MC3_07_col1_MCout[1] ^= ((cell_t)(j - 1) ^ (cell_t)j);
            colcell_tmp.k3_col = MC3_07_col1_MCout;
            colcell_tmp.SB4_B = (uint8_t) j & 0xffU;
            k3_col1_SB4_cell05[i * INIT_SB4_cell05_N + j] = colcell_tmp;
            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;

    complexity_tmp = 0ULL;
    for (size_t i = 0; i < INIT_k3_col2_N; i++)
    {
        // enumerate on cell 10 in state #MC3
        // (MC3_C[ 6], MC3_C[ 7],         *, MC3_C[ 8])
        tmp[0] = MC3_C[6];
        tmp[1] = MC3_C[7];
        tmp[2] = (uint8_t) i & 0xffU;
        tmp[3] = MC3_C[8];
        mix_column(MC3_10_col2_MCout, tmp);

        colcell_tmp.MC3_W = tmp[2];
        colcell_tmp.k3_col = MC3_10_col2_MCout;
        colcell_tmp.SB4_B = 0U;
        k3_col2_SB4_cell10[i * INIT_SB4_cell10_N] = colcell_tmp;
        for (size_t j = 1; j < INIT_SB4_cell10_N; j++)
        {
            // enumerate on cell 10 in state #SB4
            // obtain candidates for column 2 of k3 and cell 10 of #SB4
            MC3_10_col2_MCout[2] ^= ((cell_t)(j - 1) ^ (cell_t)j);
            colcell_tmp.k3_col = MC3_10_col2_MCout;
            colcell_tmp.SB4_B = (uint8_t) j & 0xffU;
            k3_col2_SB4_cell10[i * INIT_SB4_cell10_N + j] = colcell_tmp;
            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;

    complexity_tmp = 0ULL;
    for (size_t i = 0; i < INIT_k3_col3_N; i++)
    {
        // enumerate on cell 13 in state #MC3
        // (MC3_C[ 9],         *, MC3_C[10], MC3_C[11])
        tmp[0] = MC3_C[9];
        tmp[1] = (uint8_t) i & 0xffU;
        tmp[2] = MC3_C[10];
        tmp[3] = MC3_C[11];
        mix_column(MC3_13_col3_MCout, tmp);

        colcell_tmp.MC3_W = tmp[1];
        colcell_tmp.k3_col = MC3_13_col3_MCout;
        colcell_tmp.SB4_B = 0U;
        k3_col3_SB4_cell15[i * INIT_SB4_cell15_N] = colcell_tmp;
        for (size_t j = 1; j < INIT_SB4_cell15_N; j++)
        {
            // enumerate on cell 15 in state #SB4
            // obtain candidates for column 3 of k3 and cell 15 of #SB4
            MC3_13_col3_MCout[3] ^= ((cell_t)(j - 1) ^ (cell_t)j);
            colcell_tmp.k3_col = MC3_13_col3_MCout;
            colcell_tmp.SB4_B = (uint8_t) j & 0xffU;
            k3_col3_SB4_cell15[i * INIT_SB4_cell15_N + j] = colcell_tmp;
            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;

    pair<colcell_t, colcell_t > twoCols;
    column_t col0;
    column_t col1;
    column_t col2;
    column_t col3;
    colcell_t col0cell00;
    colcell_t col1cell05;
    colcell_t col2cell10;
    colcell_t col3cell15;

    cell_t k4cell05;
    cell_t k4cell10;
    cell_t k4cell15;

    complexity_tmp = 0ULL;
    for (uint32_t i = 0; i < INIT_k3_col0_SB4_cell00_N; i++)
    {
        col0 = k3_col0_SB4_cell00[i].k3_col;
        for (uint32_t j = 0; j < INIT_k3_col1_SB4_cell05_N; j++)
        {
            col1 = k3_col1_SB4_cell05[j].k3_col;
            tmp = col0 ^ col1;
            if (inv_mix_column_cell3(tmp) == MC2_C[0])
            {
                twoCols = std:: make_pair(k3_col1_SB4_cell05[j], k3_col0_SB4_cell00[i]);
                k3_col1_0.insert(twoCols);
            }
            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;
#if PRINT_COUNT
    cout << "k3_col1_0.size = 2^" << log2((double)k3_col1_0.size()) << endl;
#endif

    complexity_tmp = 0ULL;
    for (uint32_t i = 0; i < INIT_k3_col1_SB4_cell05_N; i++)
    {
        col1 = k3_col1_SB4_cell05[i].k3_col;
        for (uint32_t j = 0; j < INIT_k3_col2_SB4_cell10_N; j++)
        {
            col2 = k3_col2_SB4_cell10[j].k3_col;
            tmp = col1 ^ col2;
            if (inv_mix_column_cell2(tmp) == MC2_C[1])
            {
                twoCols = std:: make_pair(k3_col2_SB4_cell10[j], k3_col1_SB4_cell05[i]);
                k3_col2_1.insert(twoCols);
            }
            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;
#if PRINT_COUNT
    cout << "k3_col2_1.size = 2^"  << log2((double)k3_col2_1.size()) << endl;
#endif

    complexity_tmp = 0ULL;
    for (uint32_t i = 0; i < INIT_k3_col2_SB4_cell10_N; i++)
    {
        col2 = k3_col2_SB4_cell10[i].k3_col;
        for (uint32_t j = 0; j < INIT_k3_col3_SB4_cell15_N; j++)
        {
            col3 = k3_col3_SB4_cell15[j].k3_col;
            tmp = col2 ^ col3;
            if (inv_mix_column_cell1(tmp) == MC2_C[2])
            {
                twoCols = std:: make_pair(k3_col2_SB4_cell10[i], k3_col3_SB4_cell15[j]);
                k3_col2_3.insert(twoCols);
            }
            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;
#if PRINT_COUNT
    cout << "k3_col2_3.size = 2^" << log2((double)k3_col2_3.size()) << endl;
#endif

    pair<colcell_t, pair<colcell_t, colcell_t> > threeCols;
    // merge k3_col2_1 and k3_col2_3 to k3_col1_23
    complexity_tmp = 0ULL;
    for (auto& col2_1_e: k3_col2_1)
    {
        col2cell10 = col2_1_e.first;
        col1cell05 = col2_1_e.second;
        pair<MAP2_3_Iterator, MAP2_3_Iterator> col2_3_it = k3_col2_3.equal_range(col2cell10);
        for (MAP2_3_Iterator it = col2_3_it.first; it != col2_3_it.second; it++)
        {
            col3cell15 = it->second;
            twoCols = std:: make_pair(col2cell10, col3cell15);
            threeCols = std:: make_pair(col1cell05, twoCols);
            k3_col1_23.insert(threeCols);
            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;
#if PRINT_COUNT
    cout << "k3_col1_23.size = 2^"  << log2((double)k3_col1_23.size()) << endl;
#endif

    state_t SB4_B = _mm_setzero_si128();
    state_t k3_B = _mm_setzero_si128();;
    state_t k4 = _mm_setzero_si128();;
    cell_t * SB4_B_pt = (cell_t *)&SB4_B;
    cell_t * k3_B_pt = (cell_t *)&k3_B;
    cell_t * k4_pt = (cell_t *)&k4;

    uint32_t k4_C_idx = 0;
    // merge k3_col1_0 and k3_col1_23
    complexity_tmp = 0ULL;
    for (auto& col1_0_e: k3_col1_0)
    {
        col1cell05 = col1_0_e.first;
        col0cell00 = col1_0_e.second;
        pair<MAP1_23_Iterator, MAP1_23_Iterator> col1_23_it = k3_col1_23.equal_range(col1cell05);
        for (MAP1_23_Iterator it = col1_23_it.first; it != col1_23_it.second; it++)
        {
            col2cell10 = it->second.first;
            col3cell15 = it->second.second;

            k3_B = _mm_set_epi32(
                *((uint32_t *)col3cell15.k3_col.data()),
                *((uint32_t *)col2cell10.k3_col.data()),
                *((uint32_t *)col1cell05.k3_col.data()),
                *((uint32_t *)col0cell00.k3_col.data()));

            k4 = AES128_key_exp(k3_B, 4);
            k4_C_idx =   (uint32_t)k4_pt[ 5]        |
                       (((uint32_t)k4_pt[10]) << 8) |
                       (((uint32_t)k4_pt[15]) << 16);
            
            k4_C_k3SB4_B[k4_C_idx].push_back(
                ((uint64_t)(col0cell00.MC3_W)              ) |
                ((uint64_t)(col1cell05.MC3_W) << (1 * 8ULL)) |
                ((uint64_t)(col2cell10.MC3_W) << (2 * 8ULL)) |
                ((uint64_t)(col3cell15.MC3_W) << (3 * 8ULL)) |
                ((uint64_t)(col0cell00.SB4_B) << (4 * 8ULL)) |
                ((uint64_t)(col1cell05.SB4_B) << (5 * 8ULL)) |
                ((uint64_t)(col2cell10.SB4_B) << (6 * 8ULL)) |
                ((uint64_t)(col3cell15.SB4_B) << (7 * 8ULL)));

            complexity_tmp++;
        }
    }
    complexity = complexity >= complexity_tmp ? complexity : complexity_tmp;

#if PRINT_COUNT
    cout << "k4_C_k3SB4_B.size = 2^"  << log2((double)k4_C_k3SB4_B.size()) << endl;
    for (auto& a : k4_C_k3SB4_B)
    {
        cout << log2((double)a.size()) << endl;
    }
#endif
    cout << "Finish with complexity 2^" << log2((double)complexity) << endl;
}


array<vector<state_t>, AK5_CONST_SAMPLE_N> AK5_C_MC5_R;

/**************************************************************************************
 * Precompute lookup table AK5_C_MC5_R,
 * the index is the values of the 2-byte impacts on #AK5[0, 2],
 * the value is the Red cells (neutral bytes for Backward) in state #MC5
 **************************************************************************************/
void precomputeR()
{
    cout << "Start precompute the Red cells (neutral bytes for Backward) in state #MC5." << endl;
    uint64_t complexity = 0ULL;
    state_t MC5_R = _mm_setzero_si128();
    state_t AK5_R = _mm_setzero_si128();
    cell_t * MC5_R_pt = (cell_t *)&MC5_R;
    cell_t * AK5_R_pt = (cell_t *)&AK5_R;

    uint32_t AK5_C_idx = 0;

    complexity = 0ULL;
    for (uint32_t i = 0; i < (1U << 24); i++)
    {
        MC5_R_pt[1] = (cell_t)( i        & 0xffU);
        MC5_R_pt[2] = (cell_t)((i >>  8) & 0xffU);
        MC5_R_pt[3] = (cell_t)((i >> 16) & 0xffU);
        AK5_R = mix_column(MC5_R);
        AK5_C_idx = (uint32_t)AK5_R_pt[0] | ((uint32_t)AK5_R_pt[2] << 8);
        AK5_C_MC5_R[AK5_C_idx].push_back(MC5_R);
        complexity++;
    }
#if PRINT_COUNT
    for (auto& a : AK5_C_MC5_R)
    {
        cout << log2(a.size()) << endl;
    }
#endif
    cout << "Finish with complexity 2^" << log2((double)complexity) << endl;
}


/**************************************************************************************
 * Attack
 **************************************************************************************/
bool attack(state_t Target)
{
    clock_t timecnt;

    uint64_t complexity[omp_nb_threads] = {0ULL};
    uint64_t complexity_F[omp_nb_threads] = {0ULL};
    uint64_t complexity_B[omp_nb_threads] = {0ULL};
    uint64_t complexity_M[omp_nb_threads] = {0ULL};

    uint64_t complexity_Sum = 0ULL;
    uint64_t complexity_F_Sum = 0ULL;
    uint64_t complexity_B_Sum = 0ULL;
    uint64_t complexity_M_Sum = 0ULL;

    // Here, these constant impacts on MC3 and MC2 are set to be 0s.
    // One can set the impacts to be arbitrary values using block_rndfill:
    array<cell_t, 12> MC3_C = {0}; // constant impacts from Blue to Red in state #MC^3
    array<cell_t,  3> MC2_C = {0}; // constant impacts from Blue to Red in state #MC^2
    // block_rndfill(MC3_C.data(), 12);
    // block_rndfill(MC2_C.data(),  3);

    state_t Backward_C_MC3_MC = _mm_setzero_si128();
    state_t Backward_C_MC2_MC = _mm_setzero_si128();
    state_t Backward_C_MC3 = _mm_setzero_si128();
    state_t Backward_C_MC2 = _mm_setzero_si128();

    // Backward_C_MC3 and Backward_C_MC3_MC will be fixed during the whole attack
    Backward_C_MC3 = _mm_setzero_si128();
    cell_t * Backward_C_MC3_pt = (cell_t *)&Backward_C_MC3;
    Backward_C_MC3_pt[ 1] = MC3_C[ 0];
    Backward_C_MC3_pt[ 2] = MC3_C[ 1];
    Backward_C_MC3_pt[ 3] = MC3_C[ 2];
    Backward_C_MC3_pt[ 4] = MC3_C[ 3];
    Backward_C_MC3_pt[ 5] = MC3_C[ 4];
    Backward_C_MC3_pt[ 6] = MC3_C[ 5];
    Backward_C_MC3_pt[ 8] = MC3_C[ 6];
    Backward_C_MC3_pt[ 9] = MC3_C[ 7];
    Backward_C_MC3_pt[11] = MC3_C[ 8];
    Backward_C_MC3_pt[12] = MC3_C[ 9];
    Backward_C_MC3_pt[14] = MC3_C[10];
    Backward_C_MC3_pt[15] = MC3_C[11];
    Backward_C_MC3_MC = mix_column(Backward_C_MC3);

    // Backward_C_MC2 and Backward_C_MC2_MC will be fixed during the whole attack
    Backward_C_MC2 = _mm_setzero_si128();
    cell_t * Backward_C_MC2_pt = (cell_t *)&Backward_C_MC2;
    Backward_C_MC2_pt[ 7] = MC2_C[ 0];
    Backward_C_MC2_pt[10] = MC2_C[ 1];
    Backward_C_MC2_pt[13] = MC2_C[ 2];
    Backward_C_MC2_MC = mix_column(Backward_C_MC2);

    const state_t SB4_to_zero = _mm_set_epi8(State_to_Parameter_Order(
          _ , 0x52, 0x52, 0x52,
        0x52,   _ , 0x52, 0x52,
        0x52, 0x52,   _ , 0x52,
        0x52, 0x52, 0x52,   _ ));
    const state_t SB4_R_Mask = _mm_set_epi8(State_to_Parameter_Order(
          _ , 0xFF, 0xFF, 0xFF,
        0xFF,   _ , 0xFF, 0xFF,
        0xFF, 0xFF,   _ , 0xFF,
        0xFF, 0xFF, 0xFF,   _ ));
    const state_t SB5_to_zero = _mm_set_epi8(State_to_Parameter_Order(
        0x63, 0x63, 0x63, 0x63,
          _ , 0x63, 0x63, 0x63,
          _ , 0x63, 0x63, 0x63,
          _ , 0x63, 0x63, 0x63));
    const state_t MC2_C_mask = _mm_set_epi8(State_to_Parameter_Order(
          _ ,   _ ,   _ ,   _ ,
          _ ,   _ ,   _ , 0xFF,
          _ ,   _ , 0xFF,   _ ,
          _ , 0xFF,   _ ,   _ ));
    const state_t MC3_C_mask = _mm_set_epi8(State_to_Parameter_Order(
          _ , 0xFF, 0xFF, 0xFF,
        0xFF, 0xFF, 0xFF,   _ ,
        0xFF, 0xFF,   _ , 0xFF,
        0xFF,   _ , 0xFF, 0xFF));
    const state_t k4_C_mask = _mm_set_epi8(State_to_Parameter_Order(
          _ ,   _ ,   _ ,   _ ,
          _ , 0xFF,   _ ,   _ ,
          _ ,   _ , 0xFF,   _ ,
          _ ,   _ ,   _ , 0xFF));

    timecnt = clock();
    precomputeB(MC3_C, MC2_C);
    timecnt = clock() - timecnt;
    cout << "Precompute neutral bytes for forward  takes " << (((float)timecnt)/(CLOCKS_PER_SEC * 60))  << " (mins)" << endl;
    cout << endl;

    timecnt = clock();
    precomputeR();
    timecnt = clock() - timecnt;
    cout << "Precompute neutral bytes for backward takes " << (((float)timecnt)/(CLOCKS_PER_SEC * 60))  << " (mins)" << endl;
    cout << endl;

    bool findFlag = false;

    timecnt = clock();
    #pragma omp parallel for num_threads(omp_nb_threads)
    for (uint32_t k4_C_i = 0;  k4_C_i < k4_CONST_SAMPLE_N; k4_C_i++)  // sampling for the 3 Gray bytes in k4
    {
        // private variables for each thread
        int tid = omp_get_thread_num();
        clock_t timecnt_thisthread;
        state_t AK5_C = _mm_setzero_si128(); // constant impacts from Red to Blue in state #AK^5
        state_t AK4_C = _mm_setzero_si128(); // constant in state AK_4
        state_t  k4_C = _mm_setzero_si128(); // constant in state k4
        state_t AK4_k4_C_invMC = _mm_setzero_si128(); // invMixColumns(AK4_C ⊕ k4_C)
        state_t SB4_B = _mm_setzero_si128();
        state_t k3_tmp = _mm_setzero_si128();
        cell_t * k3_tmp_pt = (cell_t * )&k3_tmp;
        cell_t * AK5_C_pt = (cell_t *)&AK5_C;
        cell_t * AK4_C_pt = (cell_t *)&AK4_C;
        cell_t *  k4_C_pt = (cell_t *)& k4_C;
        cell_t * AK4_k4_C_invMC_pt = (cell_t *)& AK4_k4_C_invMC;
        cell_t * SB4_B_pt = (cell_t *)&SB4_B;

        cell_t Backward_match;
        cell_t Forward_match;
        multimap<cell_t, pair<state_t, array<state_t, 9> > > Forward_L;
        typedef std::multimap<cell_t, pair<state_t, array<state_t, 9> > >::iterator Forward_L_Iterator;
        array<state_t, 9> ks;
        vector<pair<state_t, array<state_t, 9> > > SB4_RoundKeys_B;
        vector<pair<state_t, array<state_t, 9> > >:: iterator SB4_RoundKeys_B_pt;;

        state_t Backward_R;
        cell_t * Backward_R_pt = (cell_t *)&Backward_R;
        state_t Backward_Matched;
        cell_t * Backward_Matched_pt = (cell_t *)&Backward_Matched;
        state_t zeroKey = _mm_setzero_si128();

        state_t Forward_B;
        cell_t * Forward_B_pt = (cell_t *)&Forward_B;
        state_t Forward_Matched;
        cell_t * Forward_Matched_pt = (cell_t *)&Forward_Matched;

        state_t Start_R;
        cell_t * Start_R_pt = (cell_t *)&Start_R;

        state_t kminus1;
        state_t k0;
        state_t k1_invMC;
        state_t k2;
        state_t k3;
        state_t k4;
        state_t k5;
        state_t k6;
        state_t k7;

        k4_C_pt[ 5] =  k4_C_i        & 0xffU;
        k4_C_pt[10] = (k4_C_i >>  8) & 0xffU;
        k4_C_pt[15] = (k4_C_i >> 16) & 0xffU;
        // Because Key-Schedule is independent with constants in #AK4 and impacts in #AK5,
        // given the value of the 3 Gray bytes in k4, the candidate round keys (expected to be 2^DoFF)
        // can be computed in advance and reused among loops on different constants in #AK4 and impacts in #AK5
        SB4_RoundKeys_B.clear();
        for (auto& k3SB4_B_e : k4_C_k3SB4_B[k4_C_i])
        {
            k3_tmp = _mm_setzero_si128();
            SB4_B = _mm_setzero_si128();
            k3_tmp_pt[ 0] = (cell_t)((k3SB4_B_e              ) & 0xffU);
            k3_tmp_pt[ 7] = (cell_t)((k3SB4_B_e >> (1 * 8ULL)) & 0xffU);
            k3_tmp_pt[10] = (cell_t)((k3SB4_B_e >> (2 * 8ULL)) & 0xffU);
            k3_tmp_pt[13] = (cell_t)((k3SB4_B_e >> (3 * 8ULL)) & 0xffU);
            SB4_B_pt[ 0]  = (cell_t)((k3SB4_B_e >> (4 * 8ULL)) & 0xffU);
            SB4_B_pt[ 5]  = (cell_t)((k3SB4_B_e >> (5 * 8ULL)) & 0xffU);
            SB4_B_pt[10]  = (cell_t)((k3SB4_B_e >> (6 * 8ULL)) & 0xffU);
            SB4_B_pt[15]  = (cell_t)((k3SB4_B_e >> (7 * 8ULL)) & 0xffU);
            k3_tmp = mix_column(k3_tmp);
            k3_tmp = _mm_xor_si128(k3_tmp, Backward_C_MC3_MC);
            k3_tmp = _mm_xor_si128(k3_tmp, SB4_B);

            ks[4] = k3_tmp;                       // = k3
            ks[5] = AES128_key_exp(ks[4], 4);     // = k4
            ks[6] = AES128_key_exp(ks[5], 5);     // = k5
            ks[7] = AES128_key_exp(ks[6], 6);     // = k6
            ks[8] = AES128_key_exp(ks[7], 7);     // = k7
            ks[3] = invAES128_key_exp(ks[4], 3);  // = k2
            ks[2] = invAES128_key_exp(ks[3], 2);  // = k1
            ks[1] = invAES128_key_exp(ks[2], 1);  // = k0
            ks[0] = invAES128_key_exp(ks[1], 0);  // = kminus1
            ks[2] = inv_mix_column(ks[2]);  // = k1_invMC

            #if DEBUG
            state_t k2tmp = _mm_xor_si128(ks[3], Backward_C_MC2_MC);
            k2tmp = inv_mix_column(k2tmp);
            k2tmp = _mm_and_si128(k2tmp, MC2_C_mask);
            if (!EQU(k2tmp, zeroKey)) cout << "k2 compute wrong, should has constant impact on MC2." << endl;

            state_t k4tmp = _mm_and_si128(ks[5], k4_C_mask);
            if (!EQU(k4tmp, k4_C)) cout << "k4 compute wrong, should has constant on cells 5, 10, 15." << endl;

            state_t k3tmp = _mm_xor_si128(ks[4], Backward_C_MC3_MC);
            k3tmp = _mm_xor_si128(k3tmp, SB4_B);
            k3tmp = inv_mix_column(k3tmp);
            k3tmp = _mm_and_si128(k3tmp, MC3_C_mask);
            if (!EQU(k3tmp, zeroKey)) cout << "k3 compute wrong, should has constant impact on MC3." << endl;

            state_t SB4_k3_tmp = _mm_xor_si128(SB4_B, ks[4]);
            SB4_k3_tmp = inv_mix_column(SB4_k3_tmp);
            SB4_k3_tmp = _mm_and_si128(SB4_k3_tmp, MC3_C_mask);
            if (!EQU(SB4_k3_tmp, Backward_C_MC3))
            {
                cout << "KS Wrong: inv_mix_column(SB4_B ⊕ k3_B) not equals Backward_C_MC3" << endl;
            }
            #endif
            SB4_RoundKeys_B.push_back(std:: make_pair(SB4_B, ks));
        }

        for (uint32_t AK5_C_i = 0; AK5_C_i < AK5_CONST_SAMPLE_N; AK5_C_i++)  // sampling for the 2 constant impact bytes (marked with C) in state #AK^5
        {
            AK5_C_pt[0] =  AK5_C_i       & 0xffU;
            AK5_C_pt[2] = (AK5_C_i >> 8) & 0xffU;

            for (uint32_t AK4_C_col1 = 0; AK4_C_col1 < AK4_CONST_COL_SAMPLE_N; AK4_C_col1++)  // sampling for 3 of the 9 Gray bytes in state #AK^4
            {
                AK4_C_pt[4] = (cell_t) (AK4_C_col1 & 0xffU);
                AK4_C_pt[6] = (cell_t) ((AK4_C_col1 >> 8) & 0xffU);
                AK4_C_pt[7] = (cell_t) ((AK4_C_col1 >> 16) & 0xffU);
                for (uint32_t AK4_C_col2 = 0; AK4_C_col2 < AK4_CONST_COL_SAMPLE_N; AK4_C_col2++)  // sampling for 3 of the 9 Gray bytes in state #AK^4
                {
                    AK4_C_pt[ 8] = (cell_t) (AK4_C_col2 & 0xffU);
                    AK4_C_pt[ 9] = (cell_t) ((AK4_C_col2 >> 8) & 0xffU);
                    AK4_C_pt[11] = (cell_t) ((AK4_C_col2 >> 16) & 0xffU);
                    for (uint32_t AK4_C_col3 = 0; AK4_C_col3 < AK4_CONST_COL_SAMPLE_N; AK4_C_col3++)  // sampling for 3 of the 9 Gray bytes in state #AK^4
                    {
                        AK4_C_pt[12] = (cell_t) (AK4_C_col3 & 0xffU);
                        AK4_C_pt[13] = (cell_t) ((AK4_C_col3 >> 8) & 0xffU);
                        AK4_C_pt[14] = (cell_t) ((AK4_C_col3 >> 16) & 0xffU);

                        AK4_k4_C_invMC = _mm_xor_si128(k4_C, AK4_C);
                        AK4_k4_C_invMC = inv_mix_column(AK4_k4_C_invMC);

                        // Forward compute
                        Forward_L.clear();
                        for (auto& SB4_RoundKeys_B_e: SB4_RoundKeys_B)
                        {
                            SB4_B = SB4_RoundKeys_B_e.first;
                            k3       = SB4_RoundKeys_B_e.second[4];
                            k4       = SB4_RoundKeys_B_e.second[5];
                            k5       = SB4_RoundKeys_B_e.second[6];
                            k6       = SB4_RoundKeys_B_e.second[7];
                            k7       = SB4_RoundKeys_B_e.second[8];
                            k1_invMC = SB4_RoundKeys_B_e.second[2];
                            k0       = SB4_RoundKeys_B_e.second[1];
                            kminus1  = SB4_RoundKeys_B_e.second[0];

                            Forward_B = SB4_B; // Start from #SB4 (only Blue)
                            Forward_B = _mm_or_si128(Forward_B, SB4_to_zero); // Set other cells to be invS[0] = 0x52, so that after SubBytes they turns to be 0
                            Forward_B = _mm_aesenc_si128(Forward_B, k4); // Encrypt one round to  #SB5 (before xor constant in state #AK4)
                            Forward_B = _mm_xor_si128(Forward_B, AK4_C); // #SB5 xor Constants in #AK4 (after xor constant in state #AK4)
                            Forward_B_pt[5] = Forward_B_pt[10] = Forward_B_pt[15] = 0x52; // Set other cells to be invS[0] = 0x52, so that after SubBytes they turns to be 0
                            Forward_B = _mm_aesenc_si128(Forward_B, k5); // Encrypt one round to  #SB6 (before xor constant impact from Backward)
                            Forward_B = _mm_xor_si128(Forward_B, AK5_C); // Encrypt one round to  #SB6 (after xor constant impact from Backward)
                            Forward_B = _mm_aesenc_si128(Forward_B, k6); // Encrypt one final round to  #SB7
                            Forward_B = _mm_aesenclast_si128(Forward_B, k7); // Encrypt one round to #AT before xor with Target
                            Forward_B = _mm_xor_si128(Forward_B, Target); // xor with Target get #AK_-1
                            Forward_B = _mm_xor_si128(Forward_B, kminus1); // xor with k_-1 get #SB0
                            Forward_B = _mm_aesenc_si128(Forward_B, k0); // Encrypt one round to #SB1
                            Forward_B = _mm_aesenclast_si128(Forward_B, k1_invMC); // Encrypt one final round to #MC1 ⊕ invMixColumns(k1)
                            Forward_match =  mul(0xd, Forward_B_pt[0]) ^ mul(0xe, Forward_B_pt[2]);
                            Forward_L.insert(std:: make_pair(Forward_match, SB4_RoundKeys_B_e));
                            complexity_F[tid]++;
                        }
                        #if PRINT_COUNT
                        cout << "The size of the Forward list in this loop: 2^" << log2(Forward_L.size()) << endl << endl;
                        #endif

                        uint32_t matchedcell = 0;
                        // Backward compute
                        for (auto & MC5_R: AK5_C_MC5_R[AK5_C_i])
                        {
                            Backward_R = _mm_or_si128(MC5_R, SB5_to_zero); // Set Red in MC5_R, and set other cells to be S[0], so that after invSubBytes they turns to be 0
                            Backward_R = _mm_aesdec_si128(Backward_R, zeroKey); // Decrypt one round to #MC4 (without xor constants in k4 and #AK4)
                            Backward_R = _mm_xor_si128(Backward_R, AK4_k4_C_invMC); // #MC4 xor invMixColumns(Constants in k4) xor invMixColumns(Constants in #AK4)
                            Backward_R = _mm_aesdeclast_si128(Backward_R, zeroKey); // Obtain #SB4 for computation in full match
                            Backward_R = _mm_and_si128(Backward_R, SB4_R_Mask);
                            Start_R = Backward_R;
                            Backward_R = inv_mix_column(Backward_R);  // To #MC3 without AddRoundKey k3
                            Backward_R = _mm_xor_si128(Backward_R, Backward_C_MC3);  // xor #MC3 with constant impact from Blue
                            Backward_R = _mm_aesdec_si128(Backward_R, Backward_C_MC2); // Decrypt one round to #MC2 (and xor #MC2 with constant impact from Blue)
                            Backward_R = _mm_aesdec_si128(Backward_R, zeroKey); // Decrypt one round to #MC1 (without AddRoundKey k1)
                            Backward_match =  mul(0xd, Backward_R_pt[0]) ^ mul(0xe, Backward_R_pt[2]);
                            complexity_B[tid]++;

                            pair<Forward_L_Iterator, Forward_L_Iterator> Match_it = Forward_L.equal_range(Backward_match);
                            if (Match_it.first != Match_it.second)
                            {
                                matchedcell += Forward_L.count(Backward_match);
                                for (Forward_L_Iterator it = Match_it.first; it != Match_it.second; it++)
                                {
                                    // Combine both Red and Blue in the state #SB4 to test full match on PARTIAL bits
                                    SB4_B    = it->second.first;
                                    k3       = it->second.second[4];
                                    k4       = it->second.second[5];
                                    k5       = it->second.second[6];
                                    k6       = it->second.second[7];
                                    k7       = it->second.second[8];
                                    k2       = it->second.second[3];
                                    k1_invMC = it->second.second[2];
                                    k0       = it->second.second[1];
                                    kminus1  = it->second.second[0];

                                    Backward_Matched = Start_R;
                                    Backward_Matched_pt[ 0] = SB4_B_pt[ 0];
                                    Backward_Matched_pt[ 5] = SB4_B_pt[ 5];
                                    Backward_Matched_pt[10] = SB4_B_pt[10];
                                    Backward_Matched_pt[15] = SB4_B_pt[15];
                                    Forward_Matched = Backward_Matched;

                                    Backward_Matched = _mm_xor_si128(Backward_Matched, k3);  // #SB4 ⊕ k3 => #AK3
                                    Backward_Matched = inv_mix_column(Backward_Matched); // Obtain #MC3 from #AK3
                                    Backward_Matched = _mm_aesdec_si128(Backward_Matched, inv_mix_column(k2)); // Decrypt one round to #MC2
                                    Backward_Matched = _mm_aesdec_si128(Backward_Matched, zeroKey); // Decrypt one round to #MC1 (no AddRoundkey k1)

                                    Forward_Matched = _mm_aesenc_si128(Forward_Matched, k4);  // Encrypt one round to  #SB5
                                    Forward_Matched = _mm_aesenc_si128(Forward_Matched, k5);  // Encrypt one round to  #SB6
                                    Forward_Matched = _mm_aesenc_si128(Forward_Matched, k6);  // Encrypt one round to  #SB7
                                    Forward_Matched = _mm_aesenclast_si128(Forward_Matched, k7);  // Encrypt one final round to #AT before xor with Target
                                    Forward_Matched = _mm_xor_si128(Forward_Matched, Target); // xor with Target get #AK_-1
                                    Forward_Matched = _mm_xor_si128(Forward_Matched, kminus1); // xor with k_-1 get #SB0
                                    Forward_Matched = _mm_aesenc_si128(Forward_Matched, k0); // Encrypt one round to  #SB1
                                    Forward_Matched = _mm_aesenclast_si128(Forward_Matched, k1_invMC); // Encrypt one final round to #MC1 ⊕ invMixColumns(k1)
                                    #if DEBUG
                                    cell_t Forward_match_tmp =  mul(0xd, Forward_Matched_pt[0]) ^ mul(0xe, Forward_Matched_pt[2]);
                                    if (Forward_match_tmp != Backward_match)
                                    {
                                        cout << "If run here, there must be something wrong: Forward_match_tmp != Backward_match" << endl;
                                    }
                                    #endif
                                    complexity_M[tid]++;

                                    Backward_Matched = _mm_and_si128(Backward_Matched, PARTIAL_MATCH_TARGET_MASK);
                                    Forward_Matched = _mm_and_si128(Forward_Matched, PARTIAL_MATCH_TARGET_MASK);
                                    if (EQU(Backward_Matched, Forward_Matched))
                                    {
                                        #pragma omp critical
                                        {
                                            findFlag = true;
                                            for (int ti = 0; ti < omp_nb_threads; ti++)
                                            {
                                                complexity_F_Sum += complexity_F[ti];
                                                complexity_B_Sum += complexity_B[ti];
                                                complexity_M_Sum += complexity_M[ti];
                                            }
                                            cout << "Find match on " << PARTIAL << " bits." << endl;
                                            cout << "Complexity (Total Forward ) is 2^" << log2((double)complexity_F_Sum) << endl;
                                            cout << "Complexity (Total Backward) is 2^" << log2((double)complexity_B_Sum) << endl;
                                            cout << "Complexity (After matched through MixColumns) is 2^" << log2((double)complexity_M_Sum) << endl;
                                            cout << "Complexity (Total) is 2^" << log2(((double)complexity_F_Sum + (double)complexity_B_Sum)/2.0 + (double)complexity_M_Sum) << endl;
                                            timecnt_thisthread = clock() - timecnt;
                                            cout << "Find match on " << PARTIAL << " bits takes " << (((float)timecnt_thisthread)/(CLOCKS_PER_SEC * 60))  << " (mins)" << endl;
                                            cout << endl;

                                            Backward_Matched = Start_R;
                                            Backward_Matched_pt[ 0] = SB4_B_pt[ 0];
                                            Backward_Matched_pt[ 5] = SB4_B_pt[ 5];
                                            Backward_Matched_pt[10] = SB4_B_pt[10];
                                            Backward_Matched_pt[15] = SB4_B_pt[15];
                                            Forward_Matched = Backward_Matched;

                                            cout << "==== Backward computation ====" << endl;

                                            cout << "#SB4 = " << endl;
                                            printstate(Backward_Matched);
                                            cout << "k3 = " << endl;
                                            printstate(k3);
                                            Backward_Matched = _mm_xor_si128(Backward_Matched, k3);
                                            cout << "#AK3 = " << endl;
                                            printstate(Backward_Matched);
                                            Backward_Matched = inv_mix_column(Backward_Matched);
                                            cout << "#MC3 = " << endl;
                                            printstate(Backward_Matched);
                                            cout << "k2 = " << endl;
                                            printstate(k2);
                                            Backward_Matched = _mm_aesdec_si128(Backward_Matched, inv_mix_column(k2));
                                            cout << "#MC2 = " << endl;
                                            printstate(Backward_Matched);
                                            Backward_Matched = _mm_aesdec_si128(Backward_Matched, zeroKey);
                                            cout << "invMixColumns(#SB2) = " << endl;
                                            printstate(Backward_Matched);

                                            cout << "==== Forward computation ====" << endl;
                                            cout << "#SB4 = " << endl;
                                            printstate(Forward_Matched);
                                            cout << "k4 = " << endl;
                                            printstate(k4);
                                            Forward_Matched = _mm_aesenc_si128(Forward_Matched, k4);
                                            cout << "#SB5 = " << endl;
                                            printstate(Forward_Matched);
                                            cout << "k5 = " << endl;
                                            printstate(k5);
                                            Forward_Matched = _mm_aesenc_si128(Forward_Matched, k5);
                                            cout << "#SB6 = " << endl;
                                            printstate(Forward_Matched);
                                            cout << "k6 = " << endl;
                                            printstate(k6);
                                            Forward_Matched = _mm_aesenc_si128(Forward_Matched, k6);
                                            cout << "#SB7 = " << endl;
                                            printstate(Forward_Matched);
                                            cout << "k7 = " << endl;
                                            printstate(k7);
                                            Forward_Matched = _mm_aesenclast_si128(Forward_Matched, k7);
                                            cout << "#AT = " << endl;
                                            printstate(Forward_Matched);
                                            cout << "Target = " << endl;
                                            printstate(Target);
                                            Forward_Matched = _mm_xor_si128(Forward_Matched, Target);
                                            cout << "#AK_-1 = " << endl;
                                            printstate(Forward_Matched);
                                            cout << "k_-1 = " << endl;
                                            printstate(kminus1);
                                            Forward_Matched = _mm_xor_si128(Forward_Matched, kminus1);
                                            cout << "#SB0 = " << endl;
                                            printstate(Forward_Matched);
                                            cout << "k0 = " << endl;
                                            printstate(k0);
                                            Forward_Matched = _mm_aesenc_si128(Forward_Matched, k0);
                                            cout << "#SB1 = " << endl;
                                            printstate(Forward_Matched);
                                            cout << "invMixColumns(k1) = " << endl;
                                            printstate(k1_invMC);
                                            Forward_Matched = _mm_aesenclast_si128(Forward_Matched, k1_invMC);
                                            cout << "#MC1 ⊕ invMixColumns(k1) = " << endl;
                                            printstate(Forward_Matched);

                                            exit(0);
                                        }
                                    }
                                }
                            }
                        }
                        #if PRINT_COUNT
                        cout << "Number of matched pairs through MixColumns in this loop: 2^" << log2(matchedcell) << endl << endl;
                        #endif

                        #pragma omp critical
                        if (findFlag) exit(0);
                    }
                    #pragma omp critical
                    if (findFlag) exit(0);
                }
                #pragma omp critical
                if (findFlag) exit(0);
            }
            #pragma omp critical
            if (findFlag) exit(0);
        }
        #pragma omp critical
        if (findFlag) exit(0);
    }

    if (!findFlag)
    {
        for (int ti = 0; ti < omp_nb_threads; ti++)
        {
            complexity_F_Sum += complexity_F[ti];
            complexity_B_Sum += complexity_B[ti];
            complexity_M_Sum += complexity_M[ti];
        }
        cout << "Did not find match on " << PARTIAL << " bits." << endl;
        cout << "Complexity (Total Forward ) is 2^" << log2((double)complexity_F_Sum) << endl;
        cout << "Complexity (Total Backward) is 2^" << log2((double)complexity_B_Sum) << endl;
        cout << "Complexity (After matched through MixColumns) is 2^" << log2((double)complexity_M_Sum) << endl;
        cout << "Complexity (Total) is 2^" << log2(((double)complexity_F_Sum + (double)complexity_B_Sum)/2.0 + (double)complexity_M_Sum) << endl;
        timecnt = clock() - timecnt;
        cout << "Find match on " << PARTIAL << " bits takes " << (((float)timecnt)/(CLOCKS_PER_SEC * 60))  << " (mins)" << endl;
        cout << endl;
    }

    return findFlag;
}



int main()
{
    // Here, the value of the target is set to be 0s.
    // One can set it to be arbitrary values using block_rndfill:
    state_t T = _mm_setzero_si128();
    // block_rndfill((uint8_t *)&T, sizeof(state_t));

    attack(T);

    return 0;
}