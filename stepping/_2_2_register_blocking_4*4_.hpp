#ifndef CLOSEBLAS_2_REGISTER_BLOCKING_2_2_H
#define CLOSEBLAS_2_REGISTER_BLOCKING_2_2_H

#define A(i, j) A[(i) + (j) * LDA]
#define B(i, j) B[(i) + (j) * LDB]
#define C(i, j) C[(i) + (j) * LDC]

void
dgemm_cell_ver_2_2(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);

void
scale_ver_2_2(double *C, int M, int N, int LDC, double scalar);

void
dgemm_ver_2_2(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC) {
    if (beta != 1.0) {
        scale_ver_2_2(C, M, N, LDC, beta);
    }
    int M4 = M & -4, N4 = N & -4;
    for (int i = 0; i < M4; i += 4) {
        for (int j = 0; j < N4; j += 4) {
            double c00 = C(i, j);
            double c01 = C(i, j + 1);
            double c02 = C(i, j + 2);
            double c03 = C(i, j + 3);
            double c10 = C(i + 1, j);
            double c11 = C(i + 1, j + 1);
            double c12 = C(i + 1, j + 2);
            double c13 = C(i + 1, j + 3);
            double c20 = C(i + 2, j);
            double c21 = C(i + 2, j + 1);
            double c22 = C(i + 2, j + 2);
            double c23 = C(i + 2, j + 3);
            double c30 = C(i + 3, j);
            double c31 = C(i + 3, j + 1);
            double c32 = C(i + 3, j + 2);
            double c33 = C(i + 3, j + 3);
            for (int k = 0; k < K; k++) {
                double a0 = alpha * A(i, k);
                double a1 = alpha * A(i + 1, k);
                double a2 = alpha * A(i + 2, k);
                double a3 = alpha * A(i + 3, k);
                double b0 = B(k, j);
                double b1 = B(k, j + 1);
                double b2 = B(k, j + 2);
                double b3 = B(k, j + 3);
                c00 += a0 * b0;
                c01 += a0 * b1;
                c02 += a0 * b2;
                c03 += a0 * b3;
                c10 += a1 * b0;
                c11 += a1 * b1;
                c12 += a1 * b2;
                c13 += a1 * b3;
                c20 += a2 * b0;
                c21 += a2 * b1;
                c22 += a2 * b2;
                c23 += a2 * b3;
                c30 += a3 * b0;
                c31 += a3 * b1;
                c32 += a3 * b2;
                c33 += a3 * b3;
            }
            C(i, j) = c00;
            C(i, j + 1) = c01;
            C(i, j + 2) = c02;
            C(i, j + 3) = c03;
            C(i + 1, j) = c10;
            C(i + 1, j + 1) = c11;
            C(i + 1, j + 2) = c12;
            C(i + 1, j + 3) = c13;
            C(i + 2, j) = c20;
            C(i + 2, j + 1) = c21;
            C(i + 2, j + 2) = c22;
            C(i + 2, j + 3) = c23;
            C(i + 3, j) = c30;
            C(i + 3, j + 1) = c31;
            C(i + 3, j + 2) = c32;
            C(i + 3, j + 3) = c33;
        }
    }
    if (M4 != M) {
        dgemm_cell_ver_2_2(M - M4, N, K, alpha, A + M4, LDA, B, LDB, 1.0, &C(M4, 0), LDC);
    }
    if (N4 != N) {
        dgemm_cell_ver_2_2(M4, N - N4, K, alpha, A, LDA, &B(0, N4), LDB, 1.0, &C(0, N4), LDC);
    }
}

void
dgemm_cell_ver_2_2(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC) {
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

void
scale_ver_2_2(double *C, int M, int N, int LDC, double scalar) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            C(i, j) *= scalar;
        }
    }
}

#endif //CLOSEBLAS_2_REGISTER_BLOCKING_2_2_H
