// Example program

#include "utils.h"
#define TINYPLY_IMPLEMENTATION
#include "ply.h"

using namespace std;
using namespace Eigen;
using namespace tinyply;

int main()
{
	read_ply_file("ico6.ply");

#ifdef _DEBUG
	system("PAUSE");
#endif
}

