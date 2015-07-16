#ifndef CYLINDER_H
#define CYLINDER_H

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

//This class represents a cylinder mesh
class Cylinder : public Mesh
{
private:
public:
	void init(double baseRadius, double topRadius, double height, int horizSlices, int vertStacks, vec3 color11, vec3 color12);
	
	
};

#endif