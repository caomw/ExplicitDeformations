#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <gl/glut.h>
#include <gl/GLU.H>
#include <vector>
#include "Vertex.h"
//#include "Edge.h"
#include "Logger.h"
#include "targa.h"

#define TEXT_SIZE 256	//Maximum length of the on screen message text

//Particle System class
//It builds a grid of particles (represented by Particle objects)
//These are then held together by springs using hooke's laws.  One edge corresponds to each spring.
//Hooke's law is used to represent the springs, with both a spring component and a damping comopnent
//Implicit methods are used for the integration.  This requires a solving system but allows much larger spring constants without instability
//and potentially allows more efficient implementation
//A gravity component is also present
class ParticleSystem
{
	public:
	ParticleSystem(Vertex * vertexList, int vertexCount, int * tetraList, int tetraCount, Logger * logger);
	~ParticleSystem();
	void initVBOs();
	void sendVBOs();
	void reset();
	virtual void doUpdate(double elapsedSeconds);
	void doCollisionDetectionAndResponse(double deltaT);
	void calculateNormals();
	void doRender(double videoWriteDeltaT);
	void doRender(double videoWriteDeltaT, glm::mat4 & projMatrix, glm::mat4 & modelViewMatrix);
	//UI Methods
	void increaseEarthGravity(double amount);
	void increaseStraightRestLength(double amount);
	void toggleInfoText();
	void toggleFullAmbient();
	void setWindowDimensions(int width, int height);
	void toggleAnimation();
	void toggleRenderMode();
	void toggleImageRendering();
	void doTransform();
	void printStateReport();
	void loadSpecialState();
	void setProgramObject(GLuint programObject) {this->programObject = programObject;}
	void setEyePos(glm::vec3 & eyePos);
	void setConstants(int K, int mu);

	protected:
	double halfWidth;					//Half the width of the original grid.  Used to make the grid initially be centered.
	double halfHeight;					//Half the hight of the original grid.  Used to make the grid initially be centered.
	
	//Tetrahedron data
	int * tetraList;					//Represents the tetrehdron vertex structures.  int because it contains indexes into the vertex lists.
	int numTetra;						//Number of tetrahedron in the mesh
	Vertex * orgVertices;				//Original set of vertices (undeformed)
	Vertex * defVertices;				//Deformed set of vertices
	int numVertices;					//Number of particles in the system
	vector<int> indices;				//Tetrahedral mesh indices

	vector<Vertex> floorVertices;		//Flor mesh vertices
	vector<int> floorIndices;			//Floor mesh indices
	
	double * normals;					//Array holding all face normals (if used)
	//double * vertexNormals;				//Normals for LIGHTING
	int * tetraCounts;					//Array holding number of tetrahedra adjacent to each vertex (used for normal averaging)

	//Deformation data
	double lambda;
	double mu;

	double kd;							//Damping constant (Following Choi's name of kd)
	double earthGravityValue;			//Acceleration rate for gravity

	double restitution;					//Restitution constant for collision response (controls how much the object bounces back up)
	double ud;							//Dynamic friction constant for colision response
	double us;							//Static friction constant for collision response

	double * massMatrix;				//Mass Matrix (for 3 dimensions, each set of 3 entries should be the same - mass for one particle)
	double * zeroVector;				//Vector containing all 0's
	
	double * currentForce;				//Current force

	

	int dimensionSquared;				//DIMENSION * DIMENSION occurs so frequently that a lot of computation can be saved by storing this

public:
	bool isAnimating;					//True if particles should move; false if not
protected:
	bool showInfoText;					//True if informational messages should be rendered on screen
	bool showNormals;					//True if normals should be rendered onscreen
	bool showOrg;						//True if original, undeformed vertices should be rendered on screen
	
	int * constraintParticles;			//Array containing indexes of all constrained particles)
	int numConstraints;					//Number of constrained particles
	int rows;							//Number of rows in the grid
	int cols;							//Number of columns in the grid
	

	int renderMode;						//If one will render surfaces.  If not one, will render a wireframe representation.
	char text[TEXT_SIZE];				//Used for on screen informational messages
	int windowWidth, windowHeight;		//Window width and height used for ortho mode when displaying informational text
	
	//Video generation variables (generates a series of numbered images that can be combined into a video with a tool)
	bool renderToImage;					//If true will begin rendering to series of numbered images
	uint8_t * screenBuffer;				//Buffer to hold screen capture
	uint8_t * screenRowTemp;			//Buffer to temporarily hold one row from the screen while flipping image verticaly
	int frameNumber;					//Numbered frame for series of images to be video
	char imageFileName[100];			//Name of output image file
	double timeSinceVideoWrite;			//Number of seconds elapsed since a frame was written to an image file

	Logger * logger;					//Reference to Logger class to perform all

	GLuint vboHandle[1];	  //handle to vertex buffer object for vertices
	GLuint indexVboHandle[1]; //handle to vertex buffer object for indices
	GLuint floorVboHandle[1];	  //handle to vertex buffer object for vertices
	GLuint floorIndexVboHandle[1]; //handle to vertex buffer object for indices


	GLuint programObject;				//Program object needed for shaders (notably lighting)
	GLfloat eyePos[3];		  //Position of the eye (for the camera)

	GLfloat lightAmbient[4];  //Normal ambient light setting
	GLfloat lightDiffuse[4];  //Diffuse light setting
	GLfloat lightSpecular[4]; //Specular light setting
	GLfloat lightPosition[4]; //Light position setting
	GLfloat lightColor[4];  //Light color

	GLfloat matAmbient[4];    //Material ambient setting
	GLfloat matDiffuse[4];    //Material diffuse setting
	GLfloat matSpecular[4];   //Material specular setting
	GLfloat matShininess[1];  //Material shininess setting

	GLfloat lightFullAmbient[4];  //Special "full" ambient setting for debugging

	bool ambientMode;
	
};
