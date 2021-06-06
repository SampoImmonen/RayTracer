#include "App.h"


Helmi::App::App()
{
	m_rtimage = RImage(m_height, m_width);
	initApp();
}

Helmi::App::App(int height, int width): m_height(height), m_width(width)
{
	m_rtimage = RImage(height, width);
	initApp();
}

bool Helmi::App::initApp()
{
	std::cout << "initializing app\n";
	//setup camera default position (remove when we upgrade to dynamic camera)
	float aspect_ratio = (float)(m_width) / m_height;
	glm::vec3 cam_pos(-0.5f, 1.5f, 0.0f);
	glm::vec3 cam_dir(0, -1.5f, -5.0f);
	glm::vec3 up(0, 1.0f, 0);
	m_camera = Camera(cam_pos, cam_dir, up, 90.0f, aspect_ratio);
	
	//setup default scene (for now)
	m_rrenderer.loadScene(MODELS + std::string("cornell.obj"));
	
	glm::mat4 model(1.0f);
	model = glm::scale(model, glm::vec3(m_scale));
	m_rrenderer.transformTriangles(model);
	m_rrenderer.constructBVH(SURFACE_AREA_HEURISTIC, 8);
	return false;
}

bool Helmi::App::loadScene(const std::string& modelpath)
{
	m_rrenderer.loadScene(modelpath);
	glm::mat4 model(1.0f);
	model = glm::scale(model, glm::vec3(m_scale));
	m_rrenderer.transformTriangles(model);
	//add logic to load bvh from file if one exists (much faster for large scenes)
	m_rrenderer.constructBVH(SURFACE_AREA_HEURISTIC, 8);
}

void Helmi::App::renderRT()
{
	m_rrenderer.render(m_rtimage, m_camera);
}

void Helmi::App::saveRTImage(const std::string& filepath)
{
	m_rtimage.toPPMFile(filepath);
}
