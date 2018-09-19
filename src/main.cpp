// Example program

#include "utils.h"
#define TINYPLY_IMPLEMENTATION
#include "ply.h"

using namespace std;
using namespace Eigen;
using namespace tinyply;

inline double h_(Vec n){
	return n.x() * n.y() * n.z();
}

inline double h(Vec n, double r){
	return (h_(n) -h_(-n))/2.0 + r;
}

int main()
{
	auto mesh = read_ply_file("ico8.ply");


	for(auto& v:mesh.v){
		v = v*h(v, 0.1);	
	}


	write_ply_file("output.ply", mesh);
#ifdef _DEBUG
	system("PAUSE");
#endif
}

