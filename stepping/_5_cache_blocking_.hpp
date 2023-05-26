#ifndef CLOSEBLAS_5_CACHE_BLOCKING_H
#define CLOSEBLAS_5_CACHE_BLOCKING_H

#include "immintrin.h"

#define A(i, j) A[(i) + (j) * LDA]
#define B(i, j) B[(i) + (j) * LDB]
#define C(i, j) C[(i) + (j) * LDC]

#define intrinsic_4x4_avx2 \
    a = _mm256_mul_pd(valpha, _mm256_loadu_pd(&A(i,k)));\
    b0 = _mm256_broadcast_sd(&B(k,j));\
    b1 = _mm256_broadcast_sd(&B(k,j+1));\
    b2 = _mm256_broadcast_sd(&B(k,j+2));\
    b3 = _mm256_broadcast_sd(&B(k,j+3));\
    c0 = _mm256_fmadd_pd(a,b0,c0);\
    c1 = _mm256_fmadd_pd(a,b1,c1);\
    c2 = _mm256_fmadd_pd(a,b2,c2);\
    c3 = _mm256_fmadd_pd(a,b3,c3);\
    k++;

void
macro_ver_5(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);

void
dgemm_cell_ver_5(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
                 int LDC);

void
scale_ver_5(double *C, int M, int N, int LDC, double scalar);

void
dgemm_ver_5(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
            int LDC) {
    if (beta != 1.0) {
        scale_ver_5(C, M, N, LDC, beta);
    }
    int m, n, k;
    int mg, ng, kg;
    int M_GAP = 192, N_GAP = 192, K_GAP = 192;
    for (m = 0; m < M; m += mg) {
        if (M - m > M_GAP) {
            mg = M_GAP;
        } else {
            mg = M - m;
        }
        for (n = 0; n < N; n += ng) {
            if (N - n > N_GAP) {
                ng = N_GAP;
            } else {
                ng = N - n;
            }
            for (k = 0; k < K; k += kg) {
                if (K - k > K_GAP) {
                    kg = K_GAP;
                } else {
                    kg = K - k;
                }
                macro_ver_5(mg, ng, kg, alpha, A, LDA, B, LDB, beta, C, LDC);
            }
        }
    }
}

void
macro_ver_5(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
            int LDC) {
    int i, j, k;
    int M4 = M & -4, N4 = N & -4, K4 = K & -4;
    __m256d valpha = _mm256_set1_pd(alpha);
    __m256d a, b0, b1, b2, b3;
    for (i = 0; i < M4; ++i) {
        for (j = 0; j < N4; ++j) {
            // 初始化四个小的结果矩阵
            __m256d c0 = _mm256_setzero_pd();
            __m256d c1 = _mm256_setzero_pd();
            __m256d c2 = _mm256_setzero_pd();
            __m256d c3 = _mm256_setzero_pd();
            for (k = 0; k < K4;) {
                intrinsic_4x4_avx2
                intrinsic_4x4_avx2
                intrinsic_4x4_avx2
                intrinsic_4x4_avx2
            }
            for (k = K4; k < K;) {
                intrinsic_4x4_avx2
            }
            // 保存结果
            _mm256_storeu_pd(&C(i, j), _mm256_add_pd(c0, _mm256_loadu_pd(&C(i, j))));
            _mm256_storeu_pd(&C(i, j + 1), _mm256_add_pd(c1, _mm256_loadu_pd(&C(i, j + 1))));
            _mm256_storeu_pd(&C(i, j + 2), _mm256_add_pd(c2, _mm256_loadu_pd(&C(i, j + 2))));
            _mm256_storeu_pd(&C(i, j + 3), _mm256_add_pd(c3, _mm256_loadu_pd(&C(i, j + 3))));
        }
    }
    if (M4 != M) {
        dgemm_cell_ver_5(M - M4, N, K, alpha, A + M4, LDA, B, LDB, 1.0, &C(M4, 0), LDC);
    }
    if (N4 != N) {
        dgemm_cell_ver_5(M4, N - N4, K, alpha, A, LDA, &B(0, N4), LDB, 1.0, &C(0, N4), LDC);
    }
}

void
dgemm_cell_ver_5(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
                 int LDC) {
    if (beta != 1.0) {
        scale_ver_5(C, M, N, LDC, beta);
    }
    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < N; ++n) {
            double temp = C(m, n);
            for (int k = 0; k < K; ++k) {
                temp += alpha * A(m, k) * B(k, n);
            }
            C(m, n) = temp;
        }
    }
}

void
scale_ver_5(double *C, int M, int N, int LDC, double scalar) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            C(i, j) *= scalar;
        }
    }
}

#endif //CLOSEBLAS_5_CACHE_BLOCKING_H
