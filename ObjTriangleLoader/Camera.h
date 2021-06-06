#pragma once

#include <glm/glm.hpp>
#include "math.h"
#include "Ray.h"

namespace Helmi {

class Camera {
public:

	Camera(){}

	Camera(
		glm::vec3 lookfrom,
		glm::vec3 lookat,
		glm::vec3 vup,
		float vfov, 
		float aspect_ratio
	) {
		auto theta = glm::radians(vfov);
		auto h = tan(theta / 2);
		auto viewport_height = 2.0f * h;
		auto viewport_width = aspect_ratio * viewport_height;

		auto focal_length = 1.0f;
		
		auto w = glm::normalize(lookfrom - lookat);
		auto u = glm::normalize(glm::cross(vup, w));
		auto v = glm::cross(w, u);

		origin = lookfrom;
		horizontal = viewport_width * u;
		vertical = viewport_height * v;
		lower_left_corner = origin - horizontal / 2.0f - vertical / 2.0f - w;
	}

	Ray get_ray(float u, float v) const {
		return Ray(origin, lower_left_corner + u * horizontal + v * vertical - origin);
	}

private:
	glm::vec3 origin;
	glm::vec3 lower_left_corner;
	glm::vec3 horizontal;
	glm::vec3 vertical;
};

}