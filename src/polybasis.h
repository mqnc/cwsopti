
#include "utils.h"

inline double poly(const Vec& dir, Idx exponent){
	return pow(dir.x(), exponent) * pow(dir.y(), exponent) * pow(dir.z(), exponent);
}

class PolyBasis{
public:

	const Idx nTerms;
	const std::vector<Vec> dirs;

	MatrixXd basis;
	MatrixXV diffBasis;

	PolyBasis(Idx nTerms, const std::vector<Vec>& dirs, double eps=1e-5):
		nTerms(nTerms),
		dirs(dirs)
	{
		Idx nDirs = dirs.size();

		basis = MatrixXd(nTerms, nDirs);
		diffBasis = MatrixXV(nTerms, nDirs);

		for(Idx ip=0; ip<nTerms; ip++){
			for(Idx iDir=0; iDir<nDirs; iDir++){
				auto dir = dirs[iDir];

				Vec t1 = perp1(dir);
				Vec t2 = perp2(dir);

				Idx xp = ip*2+1;

				double b0 = poly(dir.normalized(), xp);
				double db1 = poly((dir+eps*t1).normalized(), xp) - b0;
				double db2 = poly((dir+eps*t2).normalized(), xp) - b0;

				basis(ip, iDir) = b0;
				diffBasis(ip, iDir) = db1/eps*t1 + db2/eps*t2;
			}
		}
	}

};
