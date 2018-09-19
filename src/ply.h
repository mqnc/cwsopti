
# pragma once

#include "utils.h"
#include <tinyply/tinyply.h>

Mesh read_ply_file(const std::string & filepath)
{
	Mesh result;

	try
	{
		std::ifstream ss(filepath, std::ios::binary);
		if (ss.fail()) throw std::runtime_error("failed to open " + filepath);

		tinyply::PlyFile file;
		file.parse_header(ss);

		for (auto e : file.get_elements())
		{
			std::cout << "element - " << e.name << " (" << e.size << ")" << std::endl;
			for (auto p : e.properties) std::cout << "\tproperty - " << p.name << " (" << tinyply::PropertyTable[p.propertyType].str << ")" << std::endl;
		}

		std::shared_ptr<tinyply::PlyData> vertices, faces;

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

			result.v.resize(vertices->count);
			for (uint32_t i=0; i<vertices->count; i++){
				result.v[i] << vbuffer[i].x, vbuffer[i].y, vbuffer[i].z;
				result.v[i].normalize();
			}
		}
		if (vertices->t == Type::FLOAT64) {
			struct double3 { double x, y, z; };
			std::vector<double3> vbuffer(vertices->count);
			std::memcpy(vbuffer.data(), vertices->buffer.get(), vertices->buffer.size_bytes());

			result.v.resize(vertices->count);
			for (uint32_t i=0; i<vertices->count; i++){
				result.v[i] << vbuffer[i].x, vbuffer[i].y, vbuffer[i].z;
				result.v[i].normalize();
			}
		}

		struct triangle { uint32_t v1, v2, v3; };
		std::vector<triangle> fbuffer(faces->count);
		std::memcpy(fbuffer.data(), faces->buffer.get(), faces->buffer.size_bytes());

		result.f.resize(faces->count);
		for (uint32_t i=0; i<faces->count; i++){
			result.f[i] = {fbuffer[i].v1, fbuffer[i].v2, fbuffer[i].v3};
		}
	}
	catch (const std::exception & e)
	{
		std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
	}

	return result;
}
