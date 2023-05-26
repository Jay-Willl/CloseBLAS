#ifndef CLOSEBLAS_TEMPLATE_H
#define CLOSEBLAS_TEMPLATE_H

#define A(i, j) A[(i) + (j) * LDA]
#define B(i, j) B[(i) + (j) * LDB]
#define C(i, j) C[(i) + (j) * LDC]

void
dgemm_cell(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);

void
scale(double *C, int M, int N, int LDC, double scalar);

void
dgemm(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC) {
    if (beta != 1.0) {
        scale(C, M, N, LDC, beta);
    }

}

void
dgemm_cell(int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC) {
    if (beta != 1.0) {
        scale(C, M, N, LDC, beta);
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
scale(double *C, int M, int N, int LDC, double scalar) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            C(i, j) *= scalar;
        }
    }
}

#endif //CLOSEBLAS_TEMPLATE_H
