#ifndef VF_CPU_BENCH_H
#define VF_CPU_BENCH_H

#include <vector>
#include <chrono>
#include "common/Eigen/Dense"
#include <iostream>

#define mat_size 2000

class CPU_benchmark{

public:
    double measure_compute_pre_sec(){
        Eigen::MatrixXd A(mat_size, mat_size), B(mat_size, mat_size), C(mat_size, mat_size);
        A.fill(1.);
        B.fill(1.);
        C.fill(0.);

        double total_time = 0;

		auto start_time = std::chrono::system_clock::now();
		
		int num_iter = 10;
		for(int i = 0; i < num_iter; i++)
		{
			std::cout << "iter " << i << std::endl;
			C = C + A*B;
		}
		
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::system_clock::now() - start_time);
		std::cout << "Time: " << elapsed.count() << " ms" << std::endl;
		double total = std::pow(double(mat_size),  3)*num_iter ;
		std::cout << "Num operation: " << total << std::endl;
		std::cout << "Compute power: " << double(total) / elapsed.count() * 1000 << " per sec " << std::endl;
		
		std::cout << "Squared norm of result: " << C.squaredNorm() << std::endl;
    }
};

#endif // VF_CPU_BENCH_H
