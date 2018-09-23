
#include "sh/spherical_harmonics.h"

#include "utils.h"
// for stand-alone without utils.h:
// using Idx = std::size_t;
// using Vec = Eigen::Vector3d;
// using MatrixXV = Eigen::Matrix<Vec, Eigen::Dynamic, Eigen::Dynamic>;

Vec perp1(Vec v){
	return (Vec(v.y(), v.z(), -v.x()).cross(v)).normalized();
}

Vec perp2(Vec v){
	return v.cross(perp1(v));
}

class ShBasis{
public:

	const Idx degree;
	const Idx basisSize;
	const std::vector<Vec> dirs;

	MatrixXd basis;
	MatrixXV diffBasis;

	struct lm{int l; int m;};
	vector<lm> lmPairs;

	ShBasis(Idx degree, const std::vector<Vec>& dirs, double eps=1e-5):
		degree(degree),
		basisSize((degree+1)*(degree+1)),
		dirs(dirs)
	{
		Idx nDirs = dirs.size();

		basis = MatrixXd(basisSize, nDirs);
		diffBasis = MatrixXV(basisSize, nDirs);
		lmPairs = vector<lm>(basisSize);
		
		{	
			Idx i=0;
			for(int l=0; l<=degree; l++){
				for(int m=-l; m<=l; m++){
					lmPairs[i] = {l, m};
					i++;
				}
			}
		}

		for(Idx ilm=0; ilm<basisSize; ilm++){
			for(Idx iDir=0; iDir<nDirs; iDir++){
				auto dir = dirs[iDir];

				Vec t1 = perp1(dir);
				Vec t2 = perp2(dir);

				double b0  = sh::EvalSH(lmPairs[ilm].l, lmPairs[ilm].m, dir.normalized());
				double db1 = sh::EvalSH(lmPairs[ilm].l, lmPairs[ilm].m, (dir+eps*t1).normalized());
				double db2 = sh::EvalSH(lmPairs[ilm].l, lmPairs[ilm].m, (dir+eps*t2).normalized());

				basis(ilm, iDir) = b0;
				diffBasis(ilm, iDir) = db1/eps*t1 + db2/eps*t2;
			}
		}
	}

};
