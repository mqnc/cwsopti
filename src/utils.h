
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Eigen>

using Vec = Eigen::Vector3d;
using Mat = Eigen::Matrix3d;
using Tri = std::array<uint32_t, 3>;

struct Mesh{
	std::vector<Vec> v;
	std::vector<Tri> f;
};
