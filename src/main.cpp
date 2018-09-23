
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

int main()
{

	// load sphere mesh
	auto mesh = read_ply_file("ico8.ply");
	const Idx nVerts = mesh.v.size();

	// normalize data
	for(auto& v:mesh.v){v.normalize();} // float precision -> double precision
	mesh.n = mesh.v; // valid for sphere

	// random
	default_random_engine rnd;
	uniform_real_distribution<double> dist(-1,1);

	// precompute spherical harmonics basis
	auto start = std::chrono::high_resolution_clock::now();

	ShBasis shb(30, mesh.n, ShBasis::ODD);

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << elapsed.count() << "\n";

	vector<double> coeffs(shb.nHarmonics);
	for(auto& c:coeffs){
		c = dist(rnd) * 0.01;
	}

	for(Idx i=0; i<nVerts; i++){
		mesh.v[i] = mesh.n[i]*100;
		for(Idx ih=0; ih<shb.nHarmonics; ih++){
			mesh.v[i] += mesh.n[i] * shb.basis(ih, i) * coeffs[ih];
			mesh.v[i] += shb.diffBasis(ih, i) * coeffs[ih];
		}
	}

	// find smallest curvature radius
	double rCurvMin = 1e300;
	for(auto& e:mesh.e){
		double er = edgeRadius(mesh.v[e[0]], mesh.n[e[0]], mesh.v[e[1]], mesh.n[e[1]]);
		if(er < rCurvMin){rCurvMin = er;}
	}

	// shrink surface by that radius
	for(Idx i=0; i<nVerts; i++){
		mesh.v[i] -= mesh.n[i]*rCurvMin;
	}

	write_ply_file("output.ply", mesh);
#ifdef _DEBUG
	system("PAUSE");
#endif
}

