
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <unordered_map>
#include <Eigen/Eigen>

using Vec = Eigen::Vector3d;
using Mat = Eigen::Matrix3d;
using Tri = std::array<uint32_t, 3>;

struct Mesh{
	std::vector<Vec> v;
	std::vector<Tri> f;
};

using Edge = std::pair<uint32_t, uint32_t>;

std::vector<Edge> getEdges(Mesh){

}