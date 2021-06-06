#pragma once

#include "RTTriangle.h"
#include "Ray.h"

namespace Helmi {

struct RayhitResult {
	const Helmi::RTTriangle* tri;
	float t;
	float u, v;
	glm::vec3 point;
	Ray ray;

	RayhitResult(const Helmi::RTTriangle* tri,
		float t, float u, float v,
		glm::vec3 point, const Ray& ray
	) : tri(tri), t(t), u(u), v(v), point(point), ray(ray) {}

	RayhitResult():tri(nullptr), t(), u(), v(), point(), ray(){}

	inline operator bool() { return tri != nullptr; }
};
}