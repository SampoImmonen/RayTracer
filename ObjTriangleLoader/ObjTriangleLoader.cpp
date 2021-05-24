// ObjTriangleLoader.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <algorithm> 
#include <memory>

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "Ray.h"
#include "Camera.h"
#include "Triangle.h"
#include "RayhitResult.h"
#include "Bvh.h"
#include "Timer.h"

const int width = 800;
const int height = 600;

std::vector<Triangle> loadScene(const std::string& path) {
	
	std::vector<Triangle> triangles;

	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string warn;
	std::string err;

	//std::cout << MODELS << "\n";
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, path.c_str(), MODELS);

	if (!warn.empty()) {
		std::cout << warn << std::endl;
	}

	if (!err.empty()) {
		std::cerr << err << std::endl;
	}

	if (!ret) {
		exit(1);
	}

	// Loop over shapes
	for (size_t s = 0; s < shapes.size(); s++) {
		// Loop over faces(polygon)
		size_t index_offset = 0;
		
		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
			size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
			std::array<glm::vec3, 3> vertices;
			// Loop over vertices in the face.
			for (size_t v = 0; v < fv; v++) {
				// access to vertex
				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

				tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
				tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
				tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];
				vertices[v] = glm::vec3(vx, vy, vz);

				// Check if `normal_index` is zero or positive. negative = no normal data
				/*
				if (idx.normal_index >= 0) {
					tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
					tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
					tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
				}
				*/

				// Check if `texcoord_index` is zero or positive. negative = no texcoord data
				/*
				if (idx.texcoord_index >= 0) {
					tinyobj::real_t tx = attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
					tinyobj::real_t ty = attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
				}
				*/
				// Optional: vertex colors
				// tinyobj::real_t red   = attrib.colors[3*size_t(idx.vertex_index)+0];
				// tinyobj::real_t green = attrib.colors[3*size_t(idx.vertex_index)+1];
				// tinyobj::real_t blue  = attrib.colors[3*size_t(idx.vertex_index)+2];
			}
			triangles.push_back(Triangle(vertices));
			index_offset += fv;

			// per-face material
			//auto mat = shapes[s].mesh.material_ids[f];
			//std::cout << mat;
		}
	}
	return triangles;
}

glm::vec3 ray_color(const Ray& ray) {
	glm::dvec3 unit_dir = glm::normalize(ray.dir);
	auto t = 0.5 * (unit_dir.y + 1.0);
	return (1 - t)*glm::dvec3(1.0, 1.0, 1.0) + t * glm::dvec3(0.5, 0.7, 1.0);
}


class Image {
public:

	const int height;
	const int width;

	std::vector<glm::vec3> data;

	Image(int height, int width) : height(height), width(width) {
		//data.reserve(width * height);
		data.resize(width * height, glm::vec3(0.0f));
	}

	void setColor(int i, int j, glm::vec3 color) {
		//std::cout << i << " " << j << "\n";
		data.at(j * width + i) = color;
	}

	void vec3tostream(std::ostream& out, const glm::vec3& vec) {
		out << static_cast<int>(255.999 * vec.x) << ' '
			<< static_cast<int>(255.999 * vec.y) << ' '
			<< static_cast<int>(255.999 * vec.z) << '\n';
	}

	void printppm() {
		std::cout << "P3\n" << width << ' ' << height << "\n255\n";
		for (int j = height - 1; j >= 0; --j) {
			for (int i = 0; i < width; ++i) {
				vec3tostream(std::cout, data[j * width + i]);
			}
		}
	}
	//overloading []????
};


inline int longestdim(const BoundingBox& bb) {
	//returns the index of the longest dimension of a bounding box
	glm::vec3 dims = bb.max - bb.min;
	float max = 0;
	int index;
	for (int i = 0; i < 3; ++i) {
		if (dims[i] > max) {
			max = dims[i];
			index = i;
		}
	}
	return index;
}

inline BoundingBox combineBoundingBoxes(const BoundingBox& bb_1, const BoundingBox& bb_2) {
	//create a bounding box for two bounding boxes
	glm::vec3 min(std::min(bb_1.min.x, bb_2.min.x), std::min(bb_1.min.y, bb_2.min.y), std::min(bb_1.min.z, bb_2.min.z));
	glm::vec3 max(std::max(bb_1.max.x, bb_2.max.x), std::max(bb_1.max.y, bb_2.max.y), std::max(bb_1.max.z, bb_2.max.z));
	return BoundingBox(min, max);
}

inline BoundingBox centroidBox(const BoundingBox& bb, const glm::vec3& v) {
	glm::vec3 min(std::min(bb.min.x, v.x), std::min(bb.min.y, v.y), std::min(bb.min.z, v.z));
	glm::vec3 max(std::max(bb.max.x, v.x), std::max(bb.max.y, v.y), std::max(bb.max.z, v.z));
	return BoundingBox(min, max);
}

const int MAX_TRIS_PER_LEAF = 8;

class Renderer {
	int width;
	int height;

	float viewport_height;
	float viewport_width;
	float focal_lenght = 1.0f;

	glm::vec3 origin;
	glm::vec3 horizontal;
	glm::vec3 vertical;
	glm::vec3 lower_left_corner;

	std::vector<Triangle> m_triangles;

public:

	Bvh bvh;

	Renderer() : width(800),height(600), viewport_height(2.0f), viewport_width(2.0f*(800.0f/600.0f)){}
	Renderer(int width, int height): width(width), height(height) {
		viewport_height = 2.0;
		viewport_width = ((float)width / height) * viewport_height;

		origin = glm::vec3(0.0f);
		horizontal = glm::vec3(viewport_width, 0, 0);
		vertical = glm::vec3(0, viewport_height, 0);
		lower_left_corner = origin - horizontal / 2.0f - vertical / 2.0f - glm::vec3(0, 0, focal_lenght);
	}

	void set_triangles(const std::vector<Triangle>& triangles) {
		m_triangles = triangles;
	}

	void constructNode(
		std::unique_ptr<BvhNode>& N,
		std::vector<uint32_t>& indices,
		const std::vector<Triangle>& triangles,
		uint32_t start,
		uint32_t end)
	{

		std::unique_ptr<BvhNode> B = std::make_unique<BvhNode>();
		// construct BB for triangles in node;
		// also consruct BB for triangle centroids
		N->bb = BoundingBox();
		B->bb = BoundingBox();
		for (uint32_t i = start; i < end; ++i) {
			N->bb = combineBoundingBoxes(N->bb, triangles[indices[i]].boundingbox());
			B->bb = centroidBox(B->bb, triangles[indices[i]].centroid());
		}

		//set Node parameters
		N->start = start;
		N->end = end;
	
		//find longest dim
		int dim = longestdim(B->bb);

		//if not leaf partition triangles according to spatial median
		if ((end - start) > MAX_TRIS_PER_LEAF) {
			//std::cout << (end - start) << "\n";
			uint32_t* i = nullptr;
			//partition triangle indices according spatial median
			float spatialMedian = B->bb.min[dim] + (B->bb.max[dim] - B->bb.min[dim]) / 2;
			i = std::partition(&indices[start], &indices[start] + (end - start), [&triangles, spatialMedian, dim](int i) {return triangles[i].centroid()[dim] > spatialMedian; });
			
			//create child nodes
			uint32_t newEnd = std::distance(&indices[start], i);
			N->left = std::unique_ptr<BvhNode>(new BvhNode());
			constructNode(N->left, indices, triangles, N->start, N->start + newEnd);
			N->right = std::unique_ptr<BvhNode>(new BvhNode());
			constructNode(N->right, indices, triangles, N->start+newEnd, N->end);
		}

	}

	void contructBvh(const std::vector<Triangle>& triangles) {
		//initialize indices
		m_triangles = triangles;
		uint32_t num_triangles = triangles.size();
		std::vector<uint32_t> indices(num_triangles);
		for (uint32_t i = 0; i < num_triangles; i++) {
			indices[i] = i;
		}
		//init render bvh;
		bvh = Bvh();
		bvh.rootNode = std::unique_ptr<BvhNode>(new BvhNode());
		//recursively construct BVH
		constructNode(bvh.rootNode, indices, triangles, 0, num_triangles);
		bvh.indices = indices;
	}


	void bvhTraverse(std::unique_ptr<BvhNode>& N, const Ray& ray, const glm::vec3& invD, int& imin, float& tmin, float& umin, float& vmin) {
		if (AABBintersect(N->bb, ray.orig, invD) && (N->left != nullptr)) {
			bvhTraverse(N->left, ray, invD, imin, tmin, umin, vmin);
			bvhTraverse(N->right, ray, invD, imin, tmin, umin, vmin);
		}
		if ((N->left == nullptr) && (N->right == nullptr)) {
			for (int i = N->start; i < N->end; ++i) {
				float t, u, v;
				//triangle index
				uint32_t t_index = bvh.getIndex(i);
				if (m_triangles[t_index].mtIntersect(ray, t, u, v)) {
					if (t > 0.0f && t < tmin) {
						imin = t_index;
						tmin = t;
						umin = u;
						vmin = v;
					}
				}
			}
		}

	}

	RayhitResult Bvhraycast(const Ray& ray) {
		//Bvh traversal raycasting
		//first with recursion later with local stack
		float tmin = 100000.0f, umin = 0.0f, vmin = 0.0f;
		int imin = -1;

		glm::vec3 invD = glm::vec3(1 / ray.dir[0], 1 / ray.dir[1], 1 / ray.dir[2]);
		RayhitResult result;
		bvhTraverse(bvh.rootNode, ray, invD, imin, tmin, umin, vmin);
		if (imin != -1) {
			result = RayhitResult(&m_triangles[imin], tmin, umin, vmin, ray.at(tmin), ray);
		}
		return result;
	}

	RayhitResult raycast(const Ray& ray, const std::vector<Triangle>& triangles) const {

		float tmin = 1000.0f, umin = 0.0f, vmin = 0.0f;
		int imin = -1;

		//RaycastResult castresult;

		RayhitResult result;

		for (size_t i = 0; i < triangles.size(); ++i) {
			float t, u, v; // barycentric coordinates u and v (interpolation done as u*v0+v*v1+(1-u-v)*v2)
			if (triangles[i].mtIntersect(ray, t, u, v)){
				if (t > 0.0f && t < tmin) {
					imin = i;
					tmin = t;
					umin = u;
					vmin = v;
					//std::cout << "hit";
				}
			}
		}

		if (imin != -1) {
			//castresult = RaycastResult(...);
			//std::cout << "hit";
			//color = triangles[imin].normal;
			result = RayhitResult(&triangles[imin], tmin, umin, vmin, ray.at(tmin), ray);
		}

		return result;
	}

	void vec3tostream(std::ostream& out, const glm::vec3& vec) {
		out << static_cast<int>(255.999 * vec.x) << ' '
			<< static_cast<int>(255.999 * vec.y) << ' '
			<< static_cast<int>(255.999 * vec.z) << '\n';
	}

	void render(const std::vector<Triangle>& triangles) {
		
		float aspect_ratio = (float)(width) / height;
		
		glm::vec3 cam_pos(0.0f, 0.0f,-1);
		glm::vec3 cam_dir(0, 0, -5.0f);
		glm::vec3 up(0, 1.0f, 0);

		Camera cam(cam_pos, cam_dir, up, 90.0f, aspect_ratio);

		Image im(height, width);
		//std::cout << "P3\n" << width << ' ' << height << "\n255\n";
		for (int j = height - 1; j >= 0; --j) {
			//std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
			for (int i = 0; i < width; ++i) {

				//calculate ray direction through pixel
				auto u = float(i) / (width - 1);
				auto v = float(j) / (height - 1);
				Ray ray = cam.get_ray(u, v);
				
				//return ray hit result
				//RayhitResult result = raycast(ray, triangles);
				RayhitResult result = Bvhraycast(ray);
				//TODO SHADING CODE
				//simple shading
				glm::vec3 color(0.2, 0.5, 0.4);
				if (result.tri != nullptr) {
					color = glm::vec3(result.u, result.v,1-result.u-result.v);
				}

				//vec3tostream(std::cout, color);
				im.setColor(i, j, color);
			}
		}
		//im.printppm();
	}

};

void printvec3(glm::vec3 v){
	std::cout << v.x << " " << v.y << " " << v.z << "\n";
}

void triangletoworld(std::vector<Triangle>& triangles, const glm::mat4& model) {
	for (auto& t : triangles) {
		t.model_matrix(model);
	}
}

int main()
{

	std::vector<Triangle> tt = loadScene(MODELS+std::string("testscene.obj"));
	glm::vec3 pos = glm::vec3(0, 0, -5);
	float scale = 2.0f;
	glm::mat4 model(1.0f);
	model = glm::translate(model, pos);
	model = glm::scale(model, glm::vec3(scale));
	//std::cout << tt[0].vertices[0].z;
	triangletoworld(tt, model);
	//std::cout << std::max({ pos.x, pos.y, pos.z });
	//std::cout << tt.size();
	Renderer r(800, 600);
	{
		//Timer t;
		r.contructBvh(tt);
	}

	{	
		Timer t;
		r.render(tt);
	}
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
