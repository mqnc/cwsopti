
#include "sh/spherical_harmonics.h"

#include "utils.h"
// for stand-alone without utils.h:
// using Idx = std::size_t;
// using Vec = Eigen::Vector3d;
// using MatrixXV = Eigen::Matrix<Vec, Eigen::Dynamic, Eigen::Dynamic>;

/*
Vec perp1(Vec v){
	if(fabs(v.x()) <= fabs(v.y()) && fabs(v.x()) <= fabs(v.z())){
		return v.cross(Vec(1,0,0)).normalized();
	}
	if(fabs(v.y()) <= fabs(v.z())){
		return v.cross(Vec(0,1,0)).normalized();
	}
	return v.cross(Vec(0,0,1)).normalized();
}

Vec perp2(Vec v){
	return v.cross(perp1(v)).normalized();
}
*/

class ShBasis{
public:

	const Idx degree;
	Idx nHarmonics;
	const std::vector<Vec> dirs;

	MatrixXd basis;
	MatrixXV diffBasis;

	struct lm{int l; int m;};
	vector<lm> lmPairs;

	enum symmetry{
		EVEN, // only allow even degree harmonics where h(n) == h(-n)
		ODD, // only allow odd degree harmonics where h(n) == -h(-n)
		BOTH
	};

	ShBasis(Idx degree, const std::vector<Vec>& dirs, double sym=BOTH, double eps=1e-5):
		degree(degree),
		dirs(dirs)
	{
		Idx nDirs = dirs.size();
		
		nHarmonics=0;
		for(int l=0; l<=degree; l++){
			if(sym==EVEN && l%2==1 || sym==ODD && l%2==0){continue;}
			for(int m=-l; m<=l; m++){
				lmPairs.push_back({l, m});
				nHarmonics++;
			}
		}

		basis = MatrixXd(nHarmonics, nDirs);
		diffBasis = MatrixXV(nHarmonics, nDirs);

		for(Idx ih=0; ih<nHarmonics; ih++){
			for(Idx iDir=0; iDir<nDirs; iDir++){
				auto dir = dirs[iDir];

				Vec t1 = perp1(dir);
				Vec t2 = perp2(dir);

				double b0  = sh::EvalSH(lmPairs[ih].l, lmPairs[ih].m, dir.normalized());
				double db1 = sh::EvalSH(lmPairs[ih].l, lmPairs[ih].m, (dir+eps*t1).normalized()) - b0;
				double db2 = sh::EvalSH(lmPairs[ih].l, lmPairs[ih].m, (dir+eps*t2).normalized()) - b0;

				basis(ih, iDir) = b0;
				diffBasis(ih, iDir) = db1/eps*t1 + db2/eps*t2;
			}
		}
	}

};
