#include <iostream>
#include "util.hpp"
#include "exe.hpp"

int main() {
    int serial;
    std::cout << "please print the test serial number: ";
    std::cin >> serial;
    int TEST[15] = {0, 1, 2, 3, 4,
                    5, 6, 7, 8, 9,
                    10, 11, 12, 13, 14};
    int SIZE[30] = {100, 200, 300, 400, 500,
                    600, 700, 800, 900, 1000,
                    1100, 1200, 1300, 1400, 1500,
                    1600, 1700, 1800, 1900, 2000,
                    2100, 2200, 2300, 2400, 2500,
                    2600, 2700, 2800, 2900, 3000};

    int max_size = 3000;
    int num = 5;
    int serial_num;
    double *A = nullptr;
    double *B = nullptr;
    double *C = nullptr;
    double *Cx = nullptr;

    double alpha = 3.0;
    double beta = 2.0;

    A = new double[sizeof(double) * max_size * max_size];
    B = new double[sizeof(double) * max_size * max_size];
    C = new double[sizeof(double) * max_size * max_size];
    rand_matrix(A, max_size, max_size);
    rand_matrix(B, max_size, max_size);
    rand_matrix(C, max_size, max_size);

    if (serial == -1) {
        for (int i: SIZE) {
            int m, n, k;
            m = n = k = i;
            for (int s = 0; s < serial_num; s++) {
                double start = get_second();
                for (int j = 0; j < num; ++j) {
                    test_serial(serial, m, n, k, alpha, A, m, B, k, beta, C, m);
                }
                double end = get_second();
                double result = get_gflops(start, end, num, m, n, k);
                std::cout << "serial: " << serial << ", matrix size: " << i << ", performance: " << result << "GFLOPS"
                          << std::endl;
            }
        }
        delete[]A;
        delete[]B;
        delete[]C;
        return 0;
    } else {
        for (int i: SIZE) {
            int m, n, k;
            m = n = k = i;
            double start = get_second();
            for (int j = 0; j < num; ++j) {
                test_serial(serial, m, n, k, alpha, A, m, B, k, beta, C, m);
            }
            double end = get_second();
            double result = get_gflops(start, end, num, m, n, k);
            std::cout << "serial: " << serial << ", matrix size: " << i << ", performance: " << result << "GFLOPS"
                      << std::endl;
        }
    }

    delete[]A;
    delete[]B;
    delete[]C;
    return 0;
}


