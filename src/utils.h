
#pragma once

#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unordered_set>
#include <utility>

#include <Eigen/Eigen>

using Idx = std::size_t;
using Vec = Eigen::Vector3d;
using Mat = Eigen::Matrix3d;
using Edge = std::array<Idx, 2>;
using Tri = std::array<Idx, 3>;
using MatrixXd = Eigen::MatrixXd;
using MatrixXV = Eigen::Matrix<Vec, Eigen::Dynamic, Eigen::Dynamic>;

struct EdgeHash{
	std::size_t operator()(Edge const& e) const noexcept{
		return e[0] + e[1]*256;
	}
};

inline Edge ordered(Edge e){
	return e[0] < e[1] ? e : Edge({e[1], e[0]});
}

struct Mesh{
	std::vector<Vec> v;
	std::vector<Vec> n;
	std::vector<Tri> f;
	std::vector<Edge> e;

	void findEdges(){

		Idx n = f.size()/3*2;

		e.clear();
		e.reserve(n);
		std::unordered_set<Edge, EdgeHash> registry;
		registry.rehash(n*5/4);

		for(auto& face:f){
			for(char i=0; i<3; i++){
				Edge edge = ordered(Edge({face[i], face[(i+1)%3]}));
				if(registry.count(edge) == 0){
					registry.insert(edge);
					e.push_back(edge);
				}
			}
		}
	}
};

// curve radius of an edge acc. to https://computergraphics.stackexchange.com/a/1719
inline double edgeRadius(const Vec& p1, const Vec& n1, const Vec& p2, const Vec& n2){
	return (p2-p1).squaredNorm() / (n2-n1).dot(p2-p1);
}

inline double volume(const Mesh& mesh){
	double V=0;

	for(auto& f:mesh.f){
		V += mesh.v[f[0]] .dot( mesh.v[f[1]] .cross( mesh.v[f[2]] ));
	}

	return V;
}
