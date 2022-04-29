#ifndef VF_CPU_BENCH_H
#define VF_CPU_BENCH_H

#include <vector>
#include <chrono>
#include "common/Eigen/Dense"

#define mat_size 10000

class CPU_benchmark{
    inline int index(int i, int j){
        return i * mat_size + j;
    }
public:
    double measure_compute_pre_sec(){
        Eigen::Matrix<double, mat_size, mat_size> A, B, C;
        A.fill(1.);
        B.fill(1.);
        C.fill(0.);

        double total_time = 0;


//        std::vector<double>  A(mat_size* mat_size, 1);
//        std::vector<double>  B(mat_size* mat_size, 1);
//        std::vector<double>  C(mat_size* mat_size, 0);

//        // Matrix multiplication
//        for(int i = 0; i < mat_size; i++){
//            for(int j = 0; j < mat_size; j++){
//                double c = 0;
//                for(int k = 0; k < mat_size; k++){
//                    c += A[index(i,k)] * B[index(k,j)];
//                }
//                C[index(i,j)] = c;
//            }
//        }
        // Check C


        // Block matrix multiplication
    }
};

#endif // VF_CPU_BENCH_H
