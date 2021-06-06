#pragma once

#include <glm/glm.hpp>

namespace Helmi {

class Ray
{
public:
		
	glm::vec3 orig;
	glm::vec3 dir;

	Ray(){}
	Ray(const glm::vec3& origin, const glm::vec3& direction)
		:orig(origin), dir(direction)
	{}

	glm::vec3 origin() { return orig; }
	glm::vec3 direction() { return dir; }

	glm::dvec3 at(float t) const {
		return orig + t * dir;
	}
};
}
