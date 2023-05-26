#ifndef CLOSEBLAS_3_SIMD_ARM_H
#define CLOSEBLAS_3_SIMD_ARM_H

#include <arm_neon.h>

#define A(i, j) A[(i) + (j) * LDA]
#define B(i, j) B[(i) + (j) * LDB]
#define C(i, j) C[(i) + (j) * LDC]

void
dgemm_cell_ver_3_2(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);

void
scale_ver_3_2(double *C, int M, int N, int LDC, double scalar);

void
dgemm_ver_3_2(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC) {
    if (beta != 1.0) {
        scale_ver_3_2(C, M, N, LDC, beta);
    }
    int M4 = M & -4, N4 = N & -4;
    float64x2_t valpha = vdupq_n_f64(alpha);
    for (int i = 0; i < M4; ++i) {
        for (int j = 0; j < N4; ++j) {
            // 初始化四个小的结果矩阵
            // __m256d c0 = _mm256_setzero_pd();
            // __m256d c1 = _mm256_setzero_pd();
            // __m256d c2 = _mm256_setzero_pd();
            // __m256d c3 = _mm256_setzero_pd();
            float64x2_t c0 = vdupq_n_f64(0.0);
            float64x2_t c1 = vdupq_n_f64(0.0);
            float64x2_t c2 = vdupq_n_f64(0.0);
            float64x2_t c3 = vdupq_n_f64(0.0);
            for (int k = 0; k < K; k++) {
                // 加载A, 将其与alpha相乘
                // __m256d a = _mm256_mul_pd(valpha, _mm256_loadu_pd(&A(i, k)));
                float64x2_t a = vmulq_f64(valpha, vld1q_f64(&A(i, k)));
                // 加载B并拓展
                // __m256d b0 = _mm256_broadcast_sd(&B(k, j));
                // __m256d b1 = _mm256_broadcast_sd(&B(k, j + 1));
                // __m256d b2 = _mm256_broadcast_sd(&B(k, j + 2));
                // __m256d b3 = _mm256_broadcast_sd(&B(k, j + 3));
                float64x2_t b0 = vdupq_n_f64(B(k, j));
                float64x2_t b1 = vdupq_n_f64(B(k, j + 1));
                float64x2_t b2 = vdupq_n_f64(B(k, j + 2));
                float64x2_t b3 = vdupq_n_f64(B(k, j + 3));
                // 使用fma进行乘法和累加
                // c0 = _mm256_fmadd_pd(a, b0, c0);
                // c1 = _mm256_fmadd_pd(a, b1, c1);
                // c2 = _mm256_fmadd_pd(a, b2, c2);
                // c3 = _mm256_fmadd_pd(a, b3, c3);
                c0 = vfmaq_f64(c0, a, b0);
                c1 = vfmaq_f64(c1, a, b1);
                c2 = vfmaq_f64(c2, a, b2);
                c3 = vfmaq_f64(c3, a, b3);
            }
            // 保存结果
            // _mm256_storeu_pd(&C(i, j), _mm256_add_pd(c0, _mm256_loadu_pd(&C(i, j))));
            // _mm256_storeu_pd(&C(i, j + 1), _mm256_add_pd(c1, _mm256_loadu_pd(&C(i, j + 1))));
            // _mm256_storeu_pd(&C(i, j + 2), _mm256_add_pd(c2, _mm256_loadu_pd(&C(i, j + 2))));
            // _mm256_storeu_pd(&C(i, j + 3), _mm256_add_pd(c3, _mm256_loadu_pd(&C(i, j + 3))));
            vst1q_f64(&C(i, j), vaddq_f64(c0, vld1q_f64(&C(i, j))));
            vst1q_f64(&C(i, j + 1), vaddq_f64(c1, vld1q_f64(&C(i, j + 1))));
            vst1q_f64(&C(i, j + 2), vaddq_f64(c2, vld1q_f64(&C(i, j + 2))));
            vst1q_f64(&C(i, j + 3), vaddq_f64(c3, vld1q_f64(&C(i, j + 3))));
        }
    }
    if (M4 != M) {
        dgemm_cell_ver_3_2(M - M4, N, K, alpha, A + M4, LDA, B, LDB, 1.0, &C(M4, 0), LDC);
    }
    if (N4 != N) {
        dgemm_cell_ver_3_2(M4, N - N4, K, alpha, A, LDA, &B(0, N4), LDB, 1.0, &C(0, N4), LDC);
    }
}

void
dgemm_cell_ver_3_2(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC) {
    if (beta != 1.0) {
        scale_ver_3_2(C, M, N, LDC, beta);
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
scale_ver_3_2(double *C, int M, int N, int LDC, double scalar) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            C(i, j) *= scalar;
        }
    }
}

#endif //CLOSEBLAS_3_SIMD_ARM_H
