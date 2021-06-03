#pragma once

/*
file to include small utility functions

*/


#include <fstream>

#include <glm/glm.hpp>



//overload glm::vec3 for serialization
std::ofstream& operator<<(std::ofstream& os, const glm::vec3& v) {
	os << v.x << " " << v.y << "  " << v.z << " ";
	return os;
}

std::ifstream& operator>>(std::ifstream& is, glm::vec3& v) {
	is >> v.x >> v.y >> v.z;
	return is;
}


template <typename T>
void printline(T x) {
	std::cout << x << "\n";
}

inline glm::vec3 listToVec3(const float* l) {
	return glm::vec3(l[0], l[1], l[2]);
}
