#ifndef CLOSEBLAS_EXE_HPP
#define CLOSEBLAS_EXE_HPP

#include "./stepping/_1_baseline_.hpp"
#include "./stepping/_2_1_register_blocking_2*2_.hpp"
#include "./stepping/_2_2_register_blocking_4*4_.hpp"
//#include "./stepping/_3_1_simd_avx_.hpp"
#include "./stepping/_3_2_simd_neon_.hpp"
//#include "./stepping/_4_loop_unrolling_.hpp"
//#include "./stepping/_5_cache_blocking_.hpp"
//#include "./stepping/_6_packing_.hpp"


void
test_serial(int serial, int m, int n, int k, double alpha, double *A, int LDA, double *B, int LDB, double beta,
            double *C, int LDC) {
    switch (serial) {
        case 0:
            dgemm_ver_0(m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC);
            break;
        case 1:
            dgemm_ver_1(m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC);
            break;
        case 2:
            dgemm_ver_2_2(m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC);
            break;
        case 3:
            dgemm_ver_3_2(m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC);
            break;
        case 4:
            //dgemm_ver_4(m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC);
            break;
        case 5:
//            dgemm_ver_5(m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC);
            break;
        case 6:
//            dgemm_ver_6(m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC);
            break;
        case 7:
            break;
        case 8:
            break;
        case 9:
            break;
        default:
            break;
    }
}


#endif //CLOSEBLAS_EXE_HPP
