// ObjTriangleLoader.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define TINYOBJLOADER_IMPLEMENTATION

//#define PROFILE

#ifdef PROFILE
#define PROFILE_SCOPE(name) Timer t##__LINE__(name)
#define PROFILE_FUNC() PROFILE_SCOPE(__FUNCSIG__)
#else
#define PROFILE_SCOPE(name)
#endif


#include "Timer.h"
#include "App.h"



int main()
{
	
	PROFILE_SCOPE("main");
	Helmi::App app(600, 800);
	app.renderRT();
	app.saveRTImage("image.ppm");

}


