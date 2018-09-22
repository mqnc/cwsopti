
#include "utils.h"
#define TINYPLY_IMPLEMENTATION
#include "ply.h"
#include "sh/spherical_harmonics.h"

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
	Vec dh = diffOriginDistance(n);
	return (originDistance(n,r) - dh.dot(n))*n + dh;
}

// curve radius of an edge acc. to https://computergraphics.stackexchange.com/a/1719
inline double edgeRadius(const Vec& p1, const Vec& n1, const Vec& p2, const Vec& n2){
	return (p2-p1).squaredNorm() / (n2-n1).dot(p2-p1);
}

int main()
{

	Idx shDeg = 8; // degree of spherical harmonics
	Idx basisSize = (shDeg+1)*(shDeg+1); // number of basis functions

	auto mesh = read_ply_file("ico8.ply");
	const Idx nVerts = mesh.v.size();

	double r = 4*sqrt(6)/9;

	for(auto& v:mesh.v){v.normalize();} // float precision -> double precision
	mesh.n = mesh.v; // valid for sphere

	for(Idx i=0; i<nVerts; i++){
		auto& v = mesh.v[i];
		auto& n = mesh.n[i];

		v = contactPoint(n,r);
	}

	double coat = 1e300;
	for(auto& e:mesh.e){
		double er = edgeRadius(mesh.v[e[0]], mesh.n[e[0]], mesh.v[e[1]], mesh.n[e[1]]);
		if(er < coat){coat = er;}
	}


	MatrixXd basis(basisSize, nVerts);
	struct lm{int l; int m;};
	vector<lm> lmPairs(basisSize);
	{	
		Idx i=0;
		for(int l=0; l<=shDeg; l++){
			for(int m=-l; m<=l; m++){
				lmPairs[i] = {l, m};
				i++;
			}
		}
	}

	Idx u=0;
	for(Idx ilm=0; ilm<basisSize; ilm++){
		for(Idx iv=0; iv<nVerts; iv++){
			basis(ilm, iv) = EvalSH(lmPairs[ilm].l, lmPairs[ilm].m, mesh.n[iv]);
			u++;
		}
	}

	cout << u << "\n";

	write_ply_file("output.ply", mesh);
#ifdef _DEBUG
	system("PAUSE");
#endif
}

