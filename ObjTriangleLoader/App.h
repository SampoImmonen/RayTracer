#pragma once

#include <string>

#include "glm/glm.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include "Camera.h"
#include "RImage.h"
#include "RRenderer.h"




namespace Helmi {

	class App
	{
	public:
		App();
		App(int height, int width);
		~App() {};

		bool initApp();
		bool loadScene(const std::string& modelpath);
		void renderRT();
		void saveRTImage(const std::string& filepath);

	private:

		int m_width = 800, m_height = 600;
		//currenly only support single object in scene
		glm::vec3 m_scale = glm::vec3(15.0f);
		RRenderer m_rrenderer;
		Camera m_camera;
		RImage m_rtimage;
	};
}

