#ifndef SPHERE_H
#define SPHERE_H

#include <stdlib.h>
#include <gl/glew.h>
#include <gl/glut.h>
#include <gl/GL.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <vector>
#include "Vertex.h"
#include "Mesh.h"


using namespace std;
using namespace glm;

//This class represents a sphere mesh
class Sphere : public Mesh
{
private:
public:
	void init(double radius, int horizSlices, int vertStacks, vec3 color11, vec3 color12);
	vec4 getSphericalPoint(double radius, double theta, double phi);
	
};

#endif