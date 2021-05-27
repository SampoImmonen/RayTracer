#pragma once
#include <array>
#include <vector>

#include <glm/glm.hpp>
#include "Ray.h"



struct BoundingBox {
	//Class for axis aligned bounding box

	glm::vec3 min;
	glm::vec3 max;

	BoundingBox():min(glm::vec3(FLT_MAX)), max(glm::vec3(-FLT_MAX)){}
	BoundingBox(const glm::vec3& min, const glm::vec3& max) : min(min), max(max) {}

	inline float area() {
		glm::vec3 d = max - min;
		return 2.0f * (d.x * d.y + d.x * d.z + d.y * d.z);
	}

	inline glm::vec3 centroid() {
		return 0.5f * max + 0.5f * min;
	}

};


class Triangle {
public:
	std::array<glm::vec3, 3> vertices;
	glm::vec3 normal;
	glm::vec3 color = glm::vec3(1.0f, 1.0f, 0.2f);
	glm::mat3 M;
	glm::vec3 N;

	//TODO Materials

	Triangle(std::array<glm::vec3, 3> vertices) : vertices(vertices) {
		normal = glm::normalize(glm::cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));
		M[0] = vertices[1] - vertices[0];
		M[1] = vertices[2] - vertices[0];
		M[2] = normal;
		M = glm::transpose(M);
		M = glm::inverse(M);

		N = -M * vertices[0];
	}

	void calculate_normal() {
		normal = glm::normalize(glm::cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));
	}

	void model_matrix(const glm::mat4& model) {
		//apply model matrix to vertices
		for (auto& v : vertices) {
			v = glm::vec3(model * glm::vec4(v, 1.0f));
		}
		calculate_normal();
		M[0] = vertices[1] - vertices[0];
		M[1] = vertices[2] - vertices[0];
		M[2] = normal;
		M = glm::transpose(M);
		M = glm::inverse(M);

		N = -M * vertices[0];
	}

	bool intersect_woop(const Ray& ray, float& t, float& u, float& v) const {
		//from advanced computer graphics course
		glm::vec3 transformed_orig = M * ray.orig + N;
		glm::vec3 transformed_dir = M * ray.dir;

		t = -transformed_orig.z / transformed_dir.z;
		u = transformed_orig.x + transformed_dir.x * t;
		v = transformed_orig.y + transformed_dir.y * t;

		return u > .0f && v > .0f && u + v < 1.0f;
	}

	bool mtIntersect(const Ray& ray, float& t, float& u, float& v) const {

		float epsilon = 0.000000001f;

		glm::vec3 v0v1 = vertices[1] - vertices[0];
		glm::vec3 v0v2 = vertices[2] - vertices[0];

		glm::vec3 pvec = glm::cross(ray.dir, v0v2);
		float det = glm::dot(v0v1, pvec);

		if (std::abs(det) < epsilon) false;

		float invDet = 1 / det;

		glm::vec3 tvec = ray.orig - vertices[0];
		u = glm::dot(tvec, pvec) * invDet;
		if (u < 0 || u > 1) return false;

		glm::vec3 qvec = glm::cross(tvec, v0v1);
		v = glm::dot(ray.dir, qvec) * invDet;
		if (v < 0 || u + v >1) return false;

		t = glm::dot(v0v2, qvec) * invDet;

		return true;
	}

	bool naiveIntersect(const Ray& ray, float& t, float& u, float& v)  const {

		float epsilon = 0.000000001f;

		glm::vec3 v0v1 = vertices[1] - vertices[0];
		glm::vec3 v0v2 = vertices[2] - vertices[0];

		glm::vec3 N = glm::cross(v0v1, v0v2);
		float area = N.length() / 2.0f;

		float NdotRayDirection = glm::dot(N, ray.dir);
		if (std::abs(NdotRayDirection) < epsilon) {
			return false;
		}

		float d = glm::dot(N, vertices[0]);

		t = (glm::dot(N, ray.orig) + d) / (NdotRayDirection);
		if (t < 0) return false;
		glm::vec3 P = ray.orig + t * ray.dir;

		glm::vec3 C;

		glm::vec3 edge0 = vertices[1] - vertices[0];
		glm::vec3 vp0 = P - vertices[0];
		C = glm::cross(edge0, vp0);
		if (glm::dot(N, C) < 0) return false;

		glm::vec3 edge1 = vertices[2] - vertices[1];
		glm::vec3 vp1 = P - vertices[1];
		C = glm::cross(edge1, vp1);
		u = (C.length() / 2) / area;
		if (glm::dot(N, C) < 0) return false;

		glm::vec3 edge2 = vertices[0] - vertices[2];
		glm::vec3 vp2 = P - vertices[2];
		C = glm::cross(edge2, vp2);
		v = (C.length() / 2) / area;
		if (glm::dot(N, C) < 0) return false;

		return true;
	}

	inline glm::vec3 max() const {
		float x = std::max({ vertices[0].x, vertices[1].x, vertices[2].x });
		float y = std::max({ vertices[0].y, vertices[1].y, vertices[2].y });
		float z = std::max({ vertices[0].z, vertices[1].z, vertices[2].z });
		return glm::vec3(x, y, z);
	}

	inline glm::vec3 min() const {
		float x = std::min({ vertices[0].x, vertices[1].x, vertices[2].x });
		float y = std::min({ vertices[0].y, vertices[1].y, vertices[2].y });
		float z = std::min({ vertices[0].z, vertices[1].z, vertices[2].z });
		return glm::vec3(x, y, z);
	}

	inline float area() const {
		return glm::cross(vertices[1] - vertices[0], vertices[2] - vertices[0]).length() * 0.5f;
	}

	inline glm::vec3 centroid() const {
		return (vertices[0] + vertices[1] + vertices[2]) * (1.0f / 3.0f);
	}

	BoundingBox boundingbox() const {
		return BoundingBox(min(), max());
	}

};


inline bool AABBintersect(const BoundingBox& bb, const glm::vec3& orig, const glm::vec3& invD) {
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