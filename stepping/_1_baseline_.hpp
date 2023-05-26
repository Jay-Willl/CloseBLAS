#ifndef CLOSEBLAS_1_BASELINE_H
#define CLOSEBLAS_1_BASELINE_H

#include </opt/homebrew/Cellar/openblas/0.3.23/include/cblas.h>

void
dgemm_ver_1(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
            int LDC) {
    if (beta != 1.0) {
        for (int m = 0; m < M; ++m) {
            for (int n = 0; n < N; ++n) {
                C[m + n * LDC] *= beta;
            }
        }
    }
    for (int m = 0; m < M; ++m) {
        for (int n = 0; n < N; ++n) {
            double temp = C[m + n * LDC];
            for (int k = 0; k < K; ++k) {
                temp += alpha * A[m + k * LDA] * B[k + n * LDB];
            }
            C[m + n * LDC] = temp;
        }
    }
}

void
dgemm_ver_0(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
            int LDC) {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha, A, M, B, K, beta, C, M);
}

#endif //CLOSEBLAS_1_BASELINE_H
