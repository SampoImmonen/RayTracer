#pragma once

#include "Triangle.h"
#include "Ray.h"

struct RayhitResult {
	const Triangle* tri;
	float t;
	float u, v;
	glm::vec3 point;
	Ray ray;

	//Partial derivatives????
	RayhitResult(const Triangle* tri,
				float t, float u, float v,
				glm::vec3 point, const Ray& ray
		): tri(tri), t(t), u(u), v(v), point(point), ray(ray){}
	RayhitResult():tri(nullptr), t(), u(), v(), point(), ray(){}

	inline operator bool() { return tri != nullptr; }
};