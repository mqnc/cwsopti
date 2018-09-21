
#define SH 6 // degree of spherical harmonics

#include "utils.h"
#define TINYPLY_IMPLEMENTATION
#include "ply.h"

using namespace std;
using namespace Eigen;
using namespace tinyply;

// uninflated, possibly uncentered offset function
inline double h_(const Vec& n){
	return n.x() * n.y() * n.z();
}
// total differential
inline Vec Dh_(const Vec& n){
	return Vec(n.y()*n.z(), n.x()*n.z(), n.x()*n.y());
}

// inflated, centered offset function (distance between origin and tangent plane for a given direction n)
inline double h(const Vec& n, double r){
	return (h_(n) -h_(-n))/2.0 + r;
}
// total differential
inline Vec Dh(const Vec& n){
	return (Dh_(n) +Dh_(-n))/2.0;
}

// contact point for tangent plane for a given direction n
inline Vec c(const Vec& n, double r){
	Vec dh = Dh(n);
	return (h(n,r) - dh.dot(n))*n + dh;
}

// curve radius of an edge acc. to https://computergraphics.stackexchange.com/a/1719
inline double edgeRadius(const Vec& p1, const Vec& n1, const Vec& p2, const Vec& n2){
	return (p2-p1).squaredNorm() / (n2-n1).dot(p2-p1);
}

int main()
{
	auto mesh = read_ply_file("ico8.ply");

	double r = 4*sqrt(6)/9 + 10;

	for(auto& v:mesh.v){v.normalize();} // float precision -> double precision
	mesh.n = mesh.v; // valid for sphere

	for(Idx i=0; i<mesh.v.size(); i++){
		auto& v = mesh.v[i];
		auto& n = mesh.n[i];

		v = c(n,r);
	}

	double coat = 1e300;
	for(auto& e:mesh.e){
		double er = edgeRadius(mesh.v[e[0]], mesh.n[e[0]], mesh.v[e[1]], mesh.n[e[1]]);
		if(er < coat){coat = er;}
	}

	cout << coat << "\n";

	write_ply_file("output.ply", mesh);
#ifdef _DEBUG
	system("PAUSE");
#endif
}

