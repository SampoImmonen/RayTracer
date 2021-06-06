#include "RImage.h"

Helmi::RImage::RImage(int height, int width): m_height(height), m_width(width)
{
	m_data.resize(m_width * m_height, glm::vec3(0.0f));
}

void Helmi::RImage::setColor(int i, int j, const glm::vec3& color)
{
	m_data.at(j * m_width + i) = color;
}

void Helmi::RImage::toPPMFile(const std::string& filepath)
{
	std::ofstream imagefile(filepath);
	imagefile << "P3\n" << m_width << ' ' << m_height << "\n255\n";
	for (int j = m_height - 1; j >= 0; --j) {
		for (int i = 0; i < m_width; ++i) {
			vec3toStream(imagefile, m_data[j * m_width + i]);
		}
	}
}

int Helmi::RImage::getHeight() const
{
	return m_height;
}

int Helmi::RImage::getWidth() const
{
	return m_width;
}

void Helmi::RImage::vec3toStream(std::ostream& ostream, const glm::vec3& vec)
{
	ostream << static_cast<int>(255.999 * vec.x) << ' '
		<< static_cast<int>(255.999 * vec.y) << ' '
		<< static_cast<int>(255.999 * vec.z) << '\n';
}
