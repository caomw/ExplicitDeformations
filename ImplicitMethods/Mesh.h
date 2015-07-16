#ifndef MESH_H
#define MESH_H

#include <stdlib.h>
#include <gl/glew.h>
#include <gl/glut.h>
#include <gl/GL.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <vector>
#include "Vertex.h"


using namespace std;
using namespace glm;

//This class represents a mesh
//It stores the vertices and indices, and it can send it to the graphics card
class Mesh
{
protected:
	vector<VertexB> vertices;  //Mesh vertices - Use a vector so that we can take advantage of table doubling in STL
	vector<int> indices;	  //Mesh indices
	GLuint vboHandle[1];	  //handle to vertex buffer object for vertices
	GLuint indexVboHandle[1]; //handle to vertex buffer object for indices

	GLfloat lightAmbient[4];  //Normal ambient light setting
	GLfloat lightDiffuse[4];  //Diffuse light setting
	GLfloat lightSpecular[4]; //Specular light setting
	GLfloat lightPosition[4]; //Light position setting
	GLfloat lightColor[4];  //Light color

	GLfloat matAmbient[4];    //Material ambient setting
	GLfloat matDiffuse[4];    //Material diffuse setting
	GLfloat matSpecular[4];   //Material specular setting
	GLfloat matShininess[1];  //Mateiral shininess setting

	GLfloat eyePos[3];		  //Position of the eye (for the camera)

	GLfloat lightFullAmbient[4];  //Special "full" ambient setting for debugging

public:
	Mesh();
	void initVBO();
	void draw(mat4 & projMatrix, mat4 & modelViewMatrix);
	void setEyePos(vec3 & eyePos);
	void setShininess (float newShininess) { matShininess[0] = newShininess; }
};

#endif