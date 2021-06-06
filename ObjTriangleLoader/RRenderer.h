#pragma once

#include <vector>
#include <algorithm>

#include "tiny_obj_loader.h"

#include "Bvh.h"
#include "RImage.h"
#include "RTTriangle.h"
#include "Camera.h"
#include "RayhitResult.h"

//modes for bvh node creation


namespace Helmi {

	enum {
		SURFACE_AREA_HEURISTIC,
		SPATIAL_MEDIAN,
	};

	class RRenderer
	{
	public:

		RRenderer();
		~RRenderer();

		void loadScene(const std::string& filepath);
		void transformTriangles(const glm::mat4& mat);
		void constructBVH(int mode, int leaf_node_max_tris);
		void constructBVHNode(BvhNode* N, std::vector<uint16_t>& indices, const std::vector<RTTriangle>& triangles, uint16_t start, uint16_t end, int& nodecount, int mode, int leaf_node_max_tris);
		RayhitResult rayTrace(const Ray& ray);
		void render(RImage& image, const Camera& cam);
		glm::vec3 headlightShading(const RayhitResult& rt);


	private:
		std::vector<RTTriangle> m_triangles;
		Bvh m_bvh;

	};

}
