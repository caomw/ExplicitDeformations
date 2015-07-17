#include "Mesh.h"
#include "globals.h"

//Mesh constructor
//Based on http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/shading_glsl.pdf
Mesh::Mesh()
{
	lightAmbient[0] = 0.05; //Some ambient but not much at all
	lightAmbient[1] = 0.05;
	lightAmbient[2] = 0.05;
	
	lightFullAmbient[0] = 0.9; //Lots of ambient - useful for debugging in wireframe mode
	lightFullAmbient[1] = 0.9;
	lightFullAmbient[2] = 0.9;
	lightFullAmbient[3] = 1;

	lightDiffuse[0] = 0.8;
	lightDiffuse[1] = 0.8;
	lightDiffuse[2] = 0.8;
	
	lightDiffuse[3] = 1;

	lightSpecular[0] = 1;
	lightSpecular[1] = 1;
	lightSpecular[2] = 1;
	lightSpecular[3] = 1;

	lightPosition[0] = 0;
	lightPosition[1] = 0;
	lightPosition[2] = 1;
	lightPosition[3] = 1;



	matAmbient[0] = 1.0;
	matAmbient[1] = 1.0;
	matAmbient[2] = 1.0;
	matAmbient[3] = 1;

	matDiffuse[0] = 1;
	matDiffuse[1] = 1;
	matDiffuse[2] = 1;
	matDiffuse[3] = 1;

	matSpecular[0] = 0.9;
	matSpecular[1] = 0.9;
	matSpecular[2] = 0.9;
	matSpecular[3] = 1;

	//If lighting is calculated in eye space, the eye positiion basically is the origin - use this for the default
	eyePos[0] = 0.0f;
	eyePos[1] = 0.0f;
	eyePos[2] = -5.0f;

	lightColor[0] = 1.0f;
	lightColor[1] = 0.1f;
	lightColor[2] = 0.1f;
	lightColor[3] = 1.0f;

	matShininess[0] = 30;

}

//This method stores the current eye position (from the camera) for this mesh
//parameter eyePos - the eye position to store
void Mesh::setEyePos(vec3 & eyePos)
{
	this -> eyePos[0] = eyePos[0];
	this -> eyePos[1] = eyePos[1];
	this -> eyePos[2] = eyePos[2];
}

//This method initializes the VBO (gets buffers ready)
void Mesh::initVBO()
{

	glGenBuffers(1, vboHandle);
	glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(VertexB) * vertices.size(), &vertices[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glGenBuffers(1, indexVboHandle);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVboHandle[0]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * indices.size(), &indices[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

}

//This method draws the mesh
//Modified from: http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/shading_glsl.pdf
void Mesh::draw(mat4 & projMatrix, mat4 & modelViewMatrix)
{
	mat4 totalMatrix = projMatrix * modelViewMatrix;

	mat4 normalMatrix = inverse(modelViewMatrix);
	normalMatrix = transpose(normalMatrix);


	glUseProgram(programObject);

	//Parameter setup
	GLuint c0 = glGetAttribLocation(programObject, "position");
	GLuint c1 = glGetAttribLocation(programObject, "normal");
	GLuint c2 = glGetAttribLocation(programObject, "color1");

	GLuint l1 = glGetUniformLocation(programObject, "lightAmbient");
	GLuint l2 = glGetUniformLocation(programObject, "lightDiffuse");
	GLuint l3 = glGetUniformLocation(programObject, "lightSpecular");
	GLuint lp = glGetUniformLocation(programObject, "lightPosition");
	GLuint ep = glGetUniformLocation(programObject, "eyePosition");

	GLuint d1 = glGetUniformLocation(programObject, "ambient_coef");
	GLuint d2 = glGetUniformLocation(programObject, "diffuse_coef");
	GLuint d3 = glGetUniformLocation(programObject, "specular_coef");
	GLuint d4 = glGetUniformLocation(programObject, "mat_shininess");

	GLuint m1 = glGetUniformLocation(programObject, "local2clip");
	GLuint m2 = glGetUniformLocation(programObject, "local2eye");
	GLuint m3 = glGetUniformLocation(programObject, "normalMatrix");

	glEnableVertexAttribArray(c0); 
	glEnableVertexAttribArray(c1);
    glEnableVertexAttribArray(c2); 

	glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVboHandle[0]);

	glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(VertexB),(char*) NULL+0); 
	glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(VertexB),(char*) NULL+16); 
    glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(VertexB),(char*) NULL+32); 

	//If in ambient mode, add extra ambience to make things really bright
	//Otherwise use normal ambience
	if (ambientMode)
	{
		glUniform4f(l1, lightAmbient[0], lightAmbient[1], lightFullAmbient[2], 1.0);
	}
	else
	{
		glUniform4f(l1, lightAmbient[0], lightAmbient[1], lightAmbient[2], 1.0);
	}
	glUniform4f(l2, lightDiffuse[0], lightDiffuse[1], lightDiffuse[2], 1.0);
	glUniform4f(l3, lightSpecular[0], lightSpecular[1], lightSpecular[2],1.0);
	glUniform4f(lp, lightPosition[0], lightPosition[1], lightPosition[2], lightPosition[3]);
	glUniform4f(ep, eyePos[0], eyePos[1], eyePos[2], 1); 

	glUniform4f(d1, matAmbient[0], matAmbient[1], matAmbient[2], 1.0);
	glUniform4f(d2, matDiffuse[0], matDiffuse[1], matDiffuse[2], 1.0);
	glUniform4f(d3, matSpecular[0], matSpecular[1], matSpecular[2],1.0);
	glUniform1f(d4, matShininess[0]);

	glUniformMatrix4fv(m1, 1, GL_FALSE, &totalMatrix[0][0]);
	glUniformMatrix4fv(m2, 1, GL_FALSE, &modelViewMatrix[0][0]);
	glUniformMatrix4fv(m3, 1, GL_FALSE, &normalMatrix[0][0]);

	//Draw
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, (char *) NULL + 0);
	
	glUseProgram(0);


	
}