#ifndef CLOSEBLAS_6_PACKING_H
#define CLOSEBLAS_6_PACKING_H

#include "immintrin.h"

#define A(i, j) A[(i) + (j) * LDA]
#define B(i, j) B[(i) + (j) * LDB]
#define C(i, j) C[(i) + (j) * LDC]

void
dgemm_cell_ver_6(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
                 int LDC);

void
scale_ver_6(double *C, int M, int N, int LDC, double scalar);

void
pack_a_ver_6(double *src, double *dst, int dim, int i, int j);

void
pack_b_ver_6(double *src, double *dst, int dim, int i, int j);

void
dgemm_ver_6(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
            int LDC) {
    if (beta != 1.0) {
        scale_ver_6(C, M, N, LDC, beta);
    }
    int m, n, k;
    int mg, ng, kg;
    int M_GAP = 192, N_GAP = 192, K_GAP = 192;
    double *A_BUF = (double *)aligned_alloc(4096, M_GAP * K_GAP * sizeof(double));
    double *B_BUF = (double *) aligned_alloc(4096, K_GAP * N_GAP * sizeof(double));
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
            // packing B

            for (k = 0; k < K; k += kg) {
                if (K - k > K_GAP) {
                    kg = K_GAP;
                } else {
                    kg = K - k;
                }
                // packing A

                macro_ver_5(mg, ng, kg, alpha, A, LDA, B, LDB, beta, C, LDC);
            }
        }
    }
}

void
dgemm_cell_ver_6(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C,
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
scale_ver_6(double *C, int M, int N, int LDC, double scalar) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            C(i, j) *= scalar;
        }
    }
}

void
pack_a_ver_6(double *src, double *dst, int dim, int i, int j) {
    // matrix A: M * K
    double *ts, *td;
    td = dst;
    int first, second;
    int sub = i;
    for (first = 0; sub > 7; first += 8; sub -= 8) {

    }

}

void
pack_b_ver_6(double *src, double *dst, int dim, int i, int j) {
    // matrix B: K * N
    double *t1;
    double *t2;
    double *t3;
    double *t4;
    double *td;
    td = dst;
    int first, second;
    for (second = 0; second < j; second += 4) {
        t1 = src + second * dim;
        t2 = t1 + dim;
        t3 = t2 + dim;
        t4 = t3 + dim;
        for (first = 0; first < i; first++) {
            *td = *t1;
            t1++;
            td++;
            *td = *t2;
            t2++;
            td++;
            *td = *t3;
            t3++;
            td++;
            *td = *t4;
            t4++;
            td++;
        }
    }
}

#endif //CLOSEBLAS_6_PACKING_H
