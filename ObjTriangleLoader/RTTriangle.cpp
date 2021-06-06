#include "RTTriangle.h"

Helmi::RTTriangle::RTTriangle(std::array<glm::vec3, 3> vertices): m_vertices(vertices)
{
	calculateNormal();
}

Helmi::RTTriangle::RTTriangle(std::array<glm::vec3, 3> vertices, std::array<glm::vec2, 3> txCoordinates): m_vertices(vertices), m_txCoordinates(txCoordinates)
{
	calculateNormal();
}

void Helmi::RTTriangle::calculateNormal()
{
	m_normal = glm::normalize(glm::cross(m_vertices[1] - m_vertices[0], m_vertices[2] - m_vertices[0]));
}

void Helmi::RTTriangle::applyTransform(const glm::mat4& mat)
{
	for (auto& v : m_vertices) {
		v = glm::vec3(mat * glm::vec4(v, 1.0f));
	}
	calculateNormal();
}

bool Helmi::RTTriangle::intersect(const Ray& ray, float& t, float& u, float& v) const
{
	// Möller–Trumbore intersection
	float epsilon = 0.000000001f;

	glm::vec3 v0v1 = m_vertices[1] - m_vertices[0];
	glm::vec3 v0v2 = m_vertices[2] - m_vertices[0];

	glm::vec3 pvec = glm::cross(ray.dir, v0v2);
	float det = glm::dot(v0v1, pvec);

	if (std::abs(det) < epsilon) false;

	float invDet = 1 / det;

	glm::vec3 tvec = ray.orig - m_vertices[0];
	u = glm::dot(tvec, pvec) * invDet;
	if (u < 0 || u > 1) return false;

	glm::vec3 qvec = glm::cross(tvec, v0v1);
	v = glm::dot(ray.dir, qvec) * invDet;
	if (v < 0 || u + v >1) return false;

	t = glm::dot(v0v2, qvec) * invDet;

	return true;
}

inline glm::vec3 Helmi::RTTriangle::max() const
{
	float x = std::max({ m_vertices[0].x, m_vertices[1].x, m_vertices[2].x });
	float y = std::max({ m_vertices[0].y, m_vertices[1].y, m_vertices[2].y });
	float z = std::max({ m_vertices[0].z, m_vertices[1].z, m_vertices[2].z });
	return glm::vec3(x, y, z);
}

inline glm::vec3 Helmi::RTTriangle::min() const
{
	float x = std::min({ m_vertices[0].x, m_vertices[1].x, m_vertices[2].x });
	float y = std::min({ m_vertices[0].y, m_vertices[1].y, m_vertices[2].y });
	float z = std::min({ m_vertices[0].z, m_vertices[1].z, m_vertices[2].z });
	return glm::vec3(x, y, z);
}

inline float Helmi::RTTriangle::area() const
{
	return glm::cross(m_vertices[1] - m_vertices[0], m_vertices[2] - m_vertices[0]).length() * 0.5f;
}

inline glm::vec3 Helmi::RTTriangle::centroid() const
{
	return (m_vertices[0] + m_vertices[1] + m_vertices[2]) * (1.0f / 3.0f);
}

Helmi::BoundingBox Helmi::RTTriangle::boundingbox() const
{
	return BoundingBox(min(), max());
}
