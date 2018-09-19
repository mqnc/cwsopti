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

struct float2 { float x, y; };
struct float3 { float x, y, z; };
struct double3 { double x, y, z; };
struct uint3 { uint32_t x, y, z; };
struct uint4 { uint32_t x, y, z, w; };

void read_ply_file(const std::string & filepath)
{
	try
	{
		std::ifstream ss(filepath, std::ios::binary);
		if (ss.fail()) throw std::runtime_error("failed to open " + filepath);

		PlyFile file;
		file.parse_header(ss);

		std::cout << "........................................................................\n";
		for (auto c : file.get_comments()) std::cout << "Comment: " << c << std::endl;
		for (auto e : file.get_elements())
		{
			std::cout << "element - " << e.name << " (" << e.size << ")" << std::endl;
			for (auto p : e.properties) std::cout << "\tproperty - " << p.name << " (" << tinyply::PropertyTable[p.propertyType].str << ")" << std::endl;
		}
		std::cout << "........................................................................\n";

		// Tinyply treats parsed data as untyped byte buffers. See below for examples.
		std::shared_ptr<PlyData> vertices, normals, faces, texcoords;

		// The header information can be used to programmatically extract properties on elements
		// known to exist in the header prior to reading the data. For brevity of this sample, properties
		// like vertex position are hard-coded:
		try { vertices = file.request_properties_from_element("vertex", { "x", "y", "z" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { normals = file.request_properties_from_element("vertex", { "nx", "ny", "nz" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		try { texcoords = file.request_properties_from_element("vertex", { "u", "v" }); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		// Providing a list size hint (the last argument) is a 2x performance improvement. If you have
		// arbitrary ply files, it is best to leave this 0.
		try { faces = file.request_properties_from_element("face", { "vertex_indices" }, 3); }
		catch (const std::exception & e) { std::cerr << "tinyply exception: " << e.what() << std::endl; }

		file.read(ss);

		if (vertices) std::cout << "\tRead " << vertices->count << " total vertices "<< std::endl;
		if (normals) std::cout << "\tRead " << normals->count << " total vertex normals " << std::endl;
		if (texcoords) std::cout << "\tRead " << texcoords->count << " total vertex texcoords " << std::endl;
		if (faces) std::cout << "\tRead " << faces->count << " total faces (triangles) " << std::endl;

		// type casting to your own native types - Option A
		{
			const size_t numVerticesBytes = vertices->buffer.size_bytes();
			std::vector<float3> verts(vertices->count);
			std::memcpy(verts.data(), vertices->buffer.get(), numVerticesBytes);
		}

		// type casting to your own native types - Option B
		{
			std::vector<float3> verts_floats;
			std::vector<double3> verts_doubles;
			if (vertices->t == tinyply::Type::FLOAT32) { /* as floats ... */ }
			if (vertices->t == tinyply::Type::FLOAT64) { /* as doubles ... */ }
		}
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

