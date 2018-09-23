
#include "utils.h"
#define TINYPLY_IMPLEMENTATION
#include "ply.h"
#include "shbasis.h"

#define float double
#include "nelder-mead-optimizer/optimizer.h"
#undef float

using namespace std;
using namespace Eigen;
using namespace tinyply;
using namespace sh;


inline double offsetSeed(const Vec& n){
	return n.x() * n.y() * n.z();
}
inline Vec diffOffsetSeed(const Vec& n){
	return Vec(n.y()*n.z(), n.x()*n.z(), n.x()*n.y());
}

inline double originDistance(const Vec& n, double r){
	return (offsetSeed(n) -offsetSeed(-n))/2.0 + r;
}
inline Vec diffOriginDistance(const Vec& n){
	return (diffOffsetSeed(n) +diffOffsetSeed(-n))/2.0;
}

inline Vec contactPoint(const Vec& n, double r){
	
	return originDistance(n,r)*n + diffOriginDistance(n);

	// if n and diffOriginDistance(n) are not guaranteed to be perpendicular, use:
	// Vec dh = diffOriginDistance(n);
	// return (originDistance(n,r) - dh.dot(n))*n + dh;
}

// curve radius of an edge acc. to https://computergraphics.stackexchange.com/a/1719
inline double edgeRadius(const Vec& p1, const Vec& n1, const Vec& p2, const Vec& n2){
	return (p2-p1).squaredNorm() / (n2-n1).dot(p2-p1);
}

int main()
{

	// load sphere mesh
	auto mesh = read_ply_file("ico6.ply");
	const Idx nVerts = mesh.v.size();

	double r = 4*sqrt(6)/9;

	// normalize data
	for(auto& v:mesh.v){v.normalize();} // float precision -> double precision
	mesh.n = mesh.v; // valid for sphere

	for(Idx i=0; i<nVerts; i++){
		auto& v = mesh.v[i];
		auto& n = mesh.n[i];

		//v = contactPoint(n,r);
	}

	// test curve radius accuracy
	double coat = 1e300;
	for(auto& e:mesh.e){
		double er = edgeRadius(mesh.v[e[0]], mesh.n[e[0]], mesh.v[e[1]], mesh.n[e[1]]);
		if(er < coat){coat = er;}
	}

	// random
	default_random_engine rnd;
	uniform_real_distribution<double> dist(-1,1);

	// precompute spherical harmonics basis
	auto start = std::chrono::high_resolution_clock::now();

	ShBasis shb(15, mesh.n);

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << elapsed.count() << "\n";

	for(Idx ilm=0; ilm<shb.basisSize; ilm++){
		double r = dist(rnd)*0.01;

		for(Idx i=0; i<nVerts; i++){
			mesh.v[i] += mesh.n[i] * shb.basis(ilm, i) * r;
		}
	}


	write_ply_file("output.ply", mesh);
#ifdef _DEBUG
	system("PAUSE");
#endif
}

