// Example program

#include "utils.h"
#define TINYPLY_IMPLEMENTATION
#include "ply.h"

using namespace std;
using namespace Eigen;
using namespace tinyply;

// uninflated, possibly uncentered offset function
inline double h_(Vec n){
	return n.x() * n.y() * n.z();
}
// total differential
inline Vec Dh_(Vec n){
	return Vec(n.y()*n.z(), n.x()*n.z(), n.x()*n.y());
}

// inflated, centered offset function (distance between origin and tangent plane for a given direction n)
inline double h(Vec n, double r){
	return (h_(n) -h_(-n))/2.0 + r;
}
// total differential
inline Vec Dh(Vec n){
	return (Dh_(n) +Dh_(-n))/2.0;
}

// contact point for tangent plane for a given direction n
inline Vec c(Vec n, double r){
	Vec dh = Dh(n);
	return (h(n,r) - dh.dot(n))*n + dh;
}

int main()
{
	auto mesh = read_ply_file("ico8.ply");

	double r = 4*sqrt(6)/9;

	for(auto& v:mesh.v){
		v.normalize(); // float precision -> double precision
		v = c(v,r);
	}

	write_ply_file("output.ply", mesh);
#ifdef _DEBUG
	system("PAUSE");
#endif
}

