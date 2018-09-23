
#include "utils.h"
#define TINYPLY_IMPLEMENTATION
#include "ply.h"
#include "shbasis.h"
#include "nelder-mead-optimizer/optimizer.h"


using namespace std;
using namespace Eigen;
using namespace tinyply;
using namespace sh;

int main()
{

	// load sphere mesh
	auto mesh = read_ply_file("ico6.ply");
	const Idx nVerts = mesh.v.size();

	// normalize data
	for(auto& v:mesh.v){v.normalize();} // float precision -> double precision
	mesh.n = mesh.v; // valid for sphere

	// random
	default_random_engine rnd;
	uniform_real_distribution<double> dist(-1,1);

	// precompute spherical harmonics basis
	auto start = std::chrono::high_resolution_clock::now();

	ShBasis shb(15, mesh.n, ShBasis::ODD);

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	cout << elapsed.count() << "\n";

	auto eval = [&](const vector<double>& c) -> double{
		double r0 = 100; // initial radius

		// morph mesh according to coeffs using spherical harmonics
		for(Idx i=0; i<nVerts; i++){
			mesh.v[i] = mesh.n[i]*r0;
			for(Idx ih=0; ih<shb.nHarmonics; ih++){
				mesh.v[i] += mesh.n[i] * shb.basis(ih, i) * c[ih]; // normal component determined by s.h.
				mesh.v[i] += shb.diffBasis(ih, i) * c[ih]; // tangent component determined by deriv. of s.h.
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

		// compute volume per surrounding cube
		double vol = volume(mesh);
		double cubevol = pow(r0-rCurvMin, 3.0);

		// ratio is the cost we want to minimize
		return vol/cubevol;
	};

	nmopti::NelderMeadOptimizer opti(shb.nHarmonics, 0);
	vector<double> coeffs(shb.nHarmonics);

	//*
	// generate starting simplex
	for(Idx i=0; i<shb.nHarmonics+1; i++){
		for(auto& c:coeffs){c = dist(rnd) * 0.01;}
		coeffs[4]+=0.8; // a little push into tetrahedron direction
		opti.step(nmopti::Vector(coeffs), -eval(coeffs)); // minus since we have uphill simplex
		cout << i << "\n";
	}

	// do optimization
	for(Idx i=0; i<3000; i++){
		coeffs = (opti.step(nmopti::Vector(coeffs), -eval(coeffs))).coords; // minus since we have uphill simplex
		cout << i << "\n";

		if(i%100 == 0){write_ply_file("output" + to_string(i) + ".ply", mesh);}

	}
	/*/
	for(auto& c:coeffs){c = dist(rnd) * 0.01;}
	coeffs[4]+=0.8; // a little push into tetrahedron direction
	eval(coeffs);
	/**/

	write_ply_file("output.ply", mesh);
#ifdef _DEBUG
	system("PAUSE");
#endif
}

