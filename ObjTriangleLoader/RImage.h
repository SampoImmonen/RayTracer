#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <glm/glm.hpp>


namespace Helmi {

	class RImage
	{
	public:
		RImage(){}
		RImage(int height, int width);
		void setColor(int i, int j, const glm::vec3& color);
		void toPPMFile(const std::string& filepath = "image.ppm");
		int getHeight() const;
		int getWidth() const;

	private:
		int m_height, m_width;
		std::vector<glm::vec3> m_data;

		void vec3toStream(std::ostream& ostream, const glm::vec3& vec);
	};


}
