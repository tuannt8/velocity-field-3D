#pragma once
#include "Eigen/Eigen"
#include <vector>

#include <functional>

namespace VF3D {
typedef Eigen::Vector3d vec3;
typedef Eigen::Vector3i vec3i;
typedef Eigen::Matrix3d mat3;

std::vector<float> load_float_binary(std::string path);
std::vector<char> load_binary(std::string path);

void svd_decomposition(mat3 C, mat3 &L, vec3 &singularValue);

bool copy_file(std::string src, std::string dest);
}
