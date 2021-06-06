#pragma once

#include <array>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <glm/glm.hpp>
#include "Ray.h"

/*
* Class RTTriangle used during ray tracing
* Currently only triangles supprted later also support different shapes
*/


namespace Helmi {


	// Material and AABB structs used for primitives
	struct BoundingBox {
		//Class for axis aligned bounding box

		glm::vec3 min;
		glm::vec3 max;

		BoundingBox() :min(glm::vec3(FLT_MAX)), max(glm::vec3(-FLT_MAX)) {}
		BoundingBox(const glm::vec3& min, const glm::vec3& max) : min(min), max(max) {}

		inline float area() {
			glm::vec3 d = max - min;
			return 2.0f * (d.x * d.y + d.x * d.z + d.y * d.z);
		}

		inline glm::vec3 centroid() {
			return 0.5f * max + 0.5f * min;
		}

	};


	struct CMaterial {

		glm::vec3 ambient;
		glm::vec3 diffuse;
		glm::vec3 specular;
		glm::vec3 emission;

		float glossiness;

		//names of textures
		std::string ambient_tex;
		std::string diffuse_tex;
		std::string specular_tex;
		std::string bump_tex;

		CMaterial() {
			ambient = glm::vec3(0.0f);
			diffuse = glm::vec3(0.75f);
			specular = glm::vec3(0.5f);
			emission = glm::vec3(0.0f);
			glossiness = 32.0f;
		}
		~CMaterial() {}
	};


	class RTTriangle
	{
	public:

		RTTriangle(std::array<glm::vec3, 3> vertices);
		RTTriangle(std::array<glm::vec3, 3> vertices, std::array<glm::vec2, 3> txCoordinates);

		void calculateNormal();
		void applyTransform(const glm::mat4& mat);
		
		//geometric functions for Shape base class
		bool intersect(const Ray& ray, float& t, float& u, float& v) const;
		inline glm::vec3 max() const;
		inline glm::vec3 min() const;
		inline float area() const;
		inline glm::vec3 centroid() const;
		Helmi::BoundingBox boundingbox() const;


		std::array<glm::vec3, 3> m_vertices;
		std::array<glm::vec2, 3> m_txCoordinates;

		glm::vec3 m_normal;
		std::shared_ptr<Helmi::CMaterial> m_material = nullptr;
	private:
		//vertex struct to combine this info????
	};

	inline bool AABBintersect(const Helmi::BoundingBox& bb, const glm::vec3& orig, const glm::vec3& invD) {
		float start, end;
		for (int i = 0; i < 3; ++i) {
			float t0 = (bb.min[i] - orig[i]) * invD[i];
			float t1 = (bb.max[i] - orig[i]) * invD[i];
			if (invD[i] < 0.0f) {
				std::swap(t0, t1);
			}
			if (i == 0) {
				start = t0;
				end = t1;
			}
			else {
				start = t0 > start ? t0 : start;
				end = t1 < end ? t1 : end;
			}
			if (start > end) {
				return false;
			}
			if (end < 0.0f) {
				return false;
			}
		}
		return true;

	}

}

