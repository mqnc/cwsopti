// Example program
#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Eigen>
#define TINYPLY_IMPLEMENTATION
#include <tinyply/tinyply.h>

using namespace std;
using namespace Eigen;
using namespace tinyply;

using Vec = Vector3d;
using Mat = Matrix3d;

void read_ply_file(const std::string & filepath)
{
	try
	{
		std::ifstream ss(filepath, std::ios::binary);
		if (ss.fail()) throw std::runtime_error("failed to open " + filepath);

		PlyFile file;
		file.parse_header(ss);

		for (auto e : file.get_elements())
		{
			std::cout << "element - " << e.name << " (" << e.size << ")" << std::endl;
			for (auto p : e.properties) std::cout << "\tproperty - " << p.name << " (" << tinyply::PropertyTable[p.propertyType].str << ")" << std::endl;
		}

		std::shared_ptr<PlyData> vertices, faces;

		try { vertices = file.request_properties_from_element("vertex", { "x", "y", "z" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { faces = file.request_properties_from_element("face", { "vertex_indices" }, 3); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		file.read(ss);

		if (vertices) std::cout << "\tRead " << vertices->count << " total vertices "<< std::endl;
		if (faces) std::cout << "\tRead " << faces->count << " total faces (triangles) " << std::endl;

		if (vertices->t == Type::FLOAT32) {
			struct float3 { float x, y, z; };
			std::vector<float3> vbuffer(vertices->count);
			std::memcpy(vbuffer.data(), vertices->buffer.get(), vertices->buffer.size_bytes());
		}
		if (vertices->t == Type::FLOAT64) {
			struct double3 { double x, y, z; };
			std::vector<double3> vbuffer(vertices->count);
			std::memcpy(vbuffer.data(), vertices->buffer.get(), vertices->buffer.size_bytes());
		}

		struct uint3 { uint32_t x, y, z; };
		std::vector<uint3> fbuffer(faces->count);
		std::memcpy(fbuffer.data(), faces->buffer.get(), faces->buffer.size_bytes());

	}
	catch (const std::exception & e)
	{
		std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
	}
}

int main()
{
	read_ply_file("ico6.ply");
}

