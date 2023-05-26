#ifndef CLOSEBLAS_2_REGISTER_BLOCKING_2_1_H
#define CLOSEBLAS_2_REGISTER_BLOCKING_2_1_H

#define A(i, j) A[(i) + (j) * LDA]
#define B(i, j) B[(i) + (j) * LDB]
#define C(i, j) C[(i) + (j) * LDC]

void
dgemm_cell_ver_2_1(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);

void
scale_ver_2_1(double *C, int M, int N, int LDC, double scalar);

void
dgemm_ver_2_1(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC) {
    if (beta != 1.0) {
        scale_ver_2_1(C, M, N, LDC, beta);
    }
    int M2 = M & -2, N2 = N & -2;
    for (int i = 0; i < M2; i += 2) {
        for (int j = 0; j < N2; j += 2) {
            double c00 = C(i, j);
            double c01 = C(i, j + 1);
            double c10 = C(i + 1, j);
            double c11 = C(i + 1, j + 1);
            for (int k = 0; k < K; k++) {
                double a0 = alpha * A(i, k);
                double a1 = alpha * A(i + 1, k);
                double b0 = B(k, j);
                double b1 = B(k, j + 1);
                c00 += a0 * b0;
                c01 += a0 * b1;
                c10 += a1 * b0;
                c11 += a1 * b1;
            }
            C(i, j) = c00;
            C(i, j + 1) = c01;
            C(i + 1, j) = c10;
            C(i + 1, j + 1) = c11;
        }
    }
    if (M2 != M) {
        dgemm_cell_ver_2_1(M - M2, N, K, alpha, A + M2, LDA, B, LDB, 1.0, &C(M2, 0), LDC);
    }
    if (N2 != N) {
        dgemm_cell_ver_2_1(M2, N - N2, K, alpha, A, LDA, &B(0, N2), LDB, 1.0, &C(0, N2), LDC);
    }
}

void
scale_ver_2_1(double *C, int M, int N, int LDC, double scalar) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            C(i, j) *= scalar;
        }
    }
}

void
dgemm_cell_ver_2_1(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC) {
    if (beta != 1.0) {
        for (int m = 0; m < M; ++m) {
            for (int n = 0; n < N; ++n) {
                C(m, n) *= beta;
            }
        }
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

#endif //CLOSEBLAS_2_REGISTER_BLOCKING_2_1_H
