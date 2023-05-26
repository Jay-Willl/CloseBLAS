#ifndef CLOSEBLAS_UTIL_HPP
#define CLOSEBLAS_UTIL_HPP

#include <sys/time.h>
#include <ranlib.h>
#include <cstdlib>

double get_gflops(double start, double end, int num, int M, int N, int K) {
    return (2. * 1e-9 * num * M * N * K) / (end - start);
}

double get_second(){
    struct timeval time{};
    gettimeofday(&time, nullptr);
    return (time.tv_sec + 1e-6 * time.tv_usec);
}

void rand_matrix(double *A, int m, int n){
    srand(time(nullptr));
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            A[i * n + j] = (double)(rand() % 100) + 0.01 * (rand() % 100);
            if (rand() % 2 == 0) {
                A[i * n + j] *= 1.0;
            }
        }
    }
}

void copy_matrix(const double *src, double *dst, int n){
    for (int i = 0; src + i && dst + i && i < n; i++){
        *(dst + i) = *(src + i);
    }
}

#endif //CLOSEBLAS_UTIL_HPP
