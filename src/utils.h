
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
using VectorXd = Eigen::VectorXd;
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

// vector perpendicular to v
Vec perp1(const Vec& v){
	if(fabs(v.x()) <= fabs(v.y()) && fabs(v.x()) <= fabs(v.z())){
		return v.cross(Vec(1,0,0)).normalized();
	}
	if(fabs(v.y()) <= fabs(v.z())){
		return v.cross(Vec(0,1,0)).normalized();
	}
	return v.cross(Vec(0,0,1)).normalized();
}

// vector perpendicular to v and perp1(v)
Vec perp2(const Vec& v){
	return v.cross(perp1(v)).normalized();
}

struct Mesh{
	std::vector<Vec> v;
	std::vector<Vec> n;
	std::vector<Tri> f;
	std::vector<Edge> e;
	std::vector<std::vector<Idx>> nb; // vertex neighbors
	std::vector<Mat> tInv; // inverted tangential space

	void link(){

		// find edges
		Idx ne = f.size()/3*2;

		e.clear();
		e.reserve(ne);

		std::unordered_set<Edge, EdgeHash> registry;
		registry.rehash(ne*5/4);

		for(auto& face:f){
			for(char i=0; i<3; i++){
				Edge edge = ordered(Edge({face[i], face[(i+1)%3]}));
				if(registry.count(edge) == 0){
					registry.insert(edge);
					e.push_back(edge);
				}
			}
		}

		// assign neighbors
		nb.clear();
		nb.resize(v.size());

		for(auto& edge:e){
			nb[edge[0]].push_back(edge[1]);
			nb[edge[1]].push_back(edge[0]);
		}

		// compute inverse tangential spaces
		tInv.clear();
		tInv.resize(v.size());

		for(Idx i=0; i<v.size(); i++){
			tInv[i] << perp1(n[i]), perp2(n[i]), n[i];
			tInv[i].transpose();
		}

	}
};


// curve radius of an edge acc. to https://computergraphics.stackexchange.com/a/1719
inline double edgeRadius(const Vec& p1, const Vec& n1, const Vec& p2, const Vec& n2){
	return (p2-p1).squaredNorm() / (n2-n1).dot(p2-p1);
}

std::array<double, 2> curvature(const Mesh& mesh, Idx iv){
	// transform neighbor vertices into local tangential space
	// fit polynomial z = Lx^2/2 + Mxy + Ny^2/2
	// principal curvatures are eigen values of
	// / L M \
	// \ M N /
	// https://en.wikipedia.org/wiki/Second_fundamental_form

	const auto& neighbors = mesh.nb[iv];
	const Idx nnb = neighbors.size();

	MatrixXd A(nnb, 3);
	VectorXd z(nnb);

	for(Idx i=0; i<nnb; i++){
		Vec lnb = mesh.tInv[iv]*mesh.v[neighbors[i]]; // local neighbor

		A(i, 0) = lnb.x()*lnb.x()/2.0;
		A(i, 1) = lnb.x()*lnb.y();
		A(i, 2) = lnb.y()*lnb.y()/2.0;

		z(i) = lnb.z();
	}

	// solve A*x = z for x using least squares according to
	// https://eigen.tuxfamily.org/dox/group__LeastSquares.html
	Vec LMN = (A.transpose() * A).ldlt().solve(A.transpose() * z);

	// eigen values according to
	double traceHalf = (LMN[0] + LMN[2])/2.0;
	double det = LMN[0]*LMN[2] - LMN[1]*LMN[1];
	double aux = sqrt(traceHalf*traceHalf-det);
	double e1 = traceHalf + aux;
	double e2 = traceHalf - aux;

	return {min(e1, e2), max(e1, e2)};
}

inline double volume(const Mesh& mesh){
	double V=0;

	for(auto& f:mesh.f){
		V += mesh.v[f[0]] .dot( mesh.v[f[1]] .cross( mesh.v[f[2]] ));
	}

	return V;
}
