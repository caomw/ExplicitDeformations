//CSE 5542 - Real Time Rendering
//Lab 3 - Part II
//Maze generation 3D program - 1st Shaders / VBO's Implementation
//Chris Jacobsen
//Maze Generation Algorithm / Data structure based on http://mazeworks.com/mazegen/mazetut/index.htm
//Spherical coordinates based on http://mathworld.wolfram.com/SphericalCoordinates.html
//Shader usage code based on Han Wei Shen's code from class at http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/SDcubeSimple.C
//VBO code based on Han Wei Shen's code from class at http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/FFPtriangleVBOindex.cpp
//Lighting Algorithms based on http://www.cse.ohio-state.edu/~hwshen/5542/Site/Slides_files/shading_glsl.pdf and the code presented by GTA Souyma Dutta in class

#include <gl/glew.h>

//////////////////////////
//Deformations code
#include <string>
#include <sstream>
#include <iostream>
#include "StanfordSystem.h"
#include "GeorgiaInstituteSystem.h"
#include "NonlinearMethodSystem.h"
#include "ViewManager.h"
#include "Keyboard.h"
#include "TetraMeshReader.h"
//////////////////////////

#include <time.h>
#include <stdlib.h>

#include <gl/glut.h>
#include <gl/GL.h>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <list>
#include <vector>
#include <stack>
#include <iostream>
#include "Vertex.h"
#include "globals.h"
#include "Sphere.h"
#include "Cylinder.h"
#include "Cube.h"

//////
//Deformations Code
//Note: the reason these were declared globally is to accommodate Glut's function calling system
extern ParticleSystem * particleSystem;	//The main particle system
extern ViewManager viewManager;			//Instance of the view manager to allow user view control
extern Keyboard * keyboard;				//Instance of the Keyboard class to process key presses
extern Logger * logger;					//Instance of Logger class to perform all logging
extern const int whichMethod = 1;			//1 for stanford method.  2 for georgia Institute Method.  3 for NonLinear Paper method.
///////


GLuint SetupGLSL(char *fileName);

using namespace std;
using namespace glm;

vec4 getSphericalPoint(double x, double y, double z, double radius, double theta, double phi, mat4 & totalMatrix);

//Rotation / Zoom variables and constants
//http://www.cse.ohio-state.edu/~hwshen/581/Exe/main_cube.cpp
int press_x, press_y; 
int release_x, release_y; 
float x_angle = 20.0; 
float y_angle = 75.0; 
float scale_size = 1; 

int xform_mode = 0; 

#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 


//Struct to represent one point
//Used both to keep track of user's moves and to track position in recursive maze generation algorithm
struct Point
{
	int x;
	int y;
	int colorCode;
};

//Simpler version of point struct
struct PointD
{
	double x;
	double y;
};

//Struct used to represent one node in the graph for the maze
//Edges are not directly represented - it is assumed there is one edge between each pair of adjacent nodes within the graph
//This structure is based on the data structure suggested at http://mazeworks.com/mazegen/mazetut/index.htm
struct MazeCell
{
	bool topWall;	//True if have a wall at the top
	bool rightWall;
	bool bottomWall;
	bool leftWall;
	bool isVisited;	//True if this node was visited in the recursive algorithm
};

PointD startPos;	//x and y coordinates of starting point in maze
PointD endPos;	//x and y coordinates of goal point in maze
PointD userPos;	//x and y coordinates of user position in maze

const int MAZE_X_SIZE = 20;	//X dimension of maze array
const int MAZE_Y_SIZE = 20;	//Y dimension of maze array
MazeCell maze[MAZE_X_SIZE][MAZE_Y_SIZE];	//Multidimensional array representing maze - this effecively is the graph for it
double angle = 90;							//Robot's rotation angle (orientation)

double legAngle = 0;						//Curent rotation angle of robot lets
bool legForward = true;						//True if robot's legs are moving in the first direction (false if second direction)
double arm1Angle = 0;						//Main arm rotation angle
bool arm1Forward = true;					//True if robot's arm is moving in the first direction
double arm2Angle = 0;						//Sub arm rotation angle
bool arm2Forward = true;					//True if sub arm is moving in the first direction
float aspectRatio = 1.0;					//Window aspect ratio
bool fillMode = true;
bool ambientMode = false;

double yHeight = 1;							//Height of quads (y axis)
double width = 1;							//Width of one cell in maze (x axis)
double height = 1;							//Height of one cell in maze (z axis, but also used along y axis sometimes)

stack<mat4> matrixStack; //Stack to allow hiearchical pushing and popping of transform matrices

//Shader program object
GLuint programObject;

//Mesh objects
Sphere standSphere;
Cylinder cylinder;
Cube blueCube;
Cube greenCube;
Cube greenCube2;

//Mouse interaction - click response
//Modified from http://www.cse.ohio-state.edu/~hwshen/581/Exe/main_cube.cpp
void mymouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    press_x = x; press_y = y; 
    if (button == GLUT_LEFT_BUTTON)
      xform_mode = XFORM_ROTATE; 
	 else if (button == GLUT_RIGHT_BUTTON) 
      xform_mode = XFORM_SCALE; 
  }
  else if (state == GLUT_UP) {
	  xform_mode = XFORM_NONE; 
  }
}

//Mouse interaction - mouse moving
//Modified from http://www.cse.ohio-state.edu/~hwshen/581/Exe/main_cube.cpp
void mymotion(int x, int y)
{
    if (xform_mode==XFORM_ROTATE) {
      x_angle -= (x - press_x)/5.0; 
	  
	  if (x_angle > 180) x_angle -= 360; 
      else if (x_angle <-180) x_angle += 360; 
      press_x = x; 
	   
      y_angle -= (y - press_y)/5.0;
	  //When rotating the eye position, it will flip if you rotate too far on the y axis
	  if (y_angle < 0.5)
	  {
		  y_angle = 0.5;
	  }
	  else if (y_angle > 179.5)
	  {
		  y_angle = 179.5;
	  }
      
      press_y = y; 
    }
	else if (xform_mode == XFORM_SCALE){
      float old_size = scale_size;
      scale_size *= (1+ (y - press_y)/60.0); 
      if (scale_size <0.5) scale_size = old_size; 
      press_y = y; 
    }
	glutPostRedisplay(); 
}


//This function resets all values within the nodes of the array so that a new maze can be generated
void initMaze()
{
	for (int i = 0; i < MAZE_Y_SIZE; i++)
	{
		for (int j = 0; j < MAZE_X_SIZE; j++)
		{
			maze[i][j].topWall = maze[i][j].rightWall = maze[i][j].bottomWall = maze[i][j].leftWall = true;
			maze[i][j].isVisited = false;
		}
	}
}

//This function generates a new maze at random
//For a new maze, initMaze() should have been called first to clean all fields in the nodes
//Also assumes that random numbers have been seeded
//Inputs: curPos - starting point for maze generation
//Based on the algorithm at http://mazeworks.com/mazegen/mazetut/index.htm
void generateMaze(Point curPos)
{
	maze[curPos.y][curPos.x].isVisited = true;

	//This vector represents the remaining choices in our randomized depth first search
	vector <int> choices(4);
	choices[0] = 1;
	choices[1] = 2;
	choices[2] = 3;
	choices[3] = 4;

	Point newPos;	//The new node being examined in the DFS
	int choice = 0;	//The random choice we make in our randomized DFS
	while (choices.size() > 0)
	{
		//Make a random choice out of the unexplored nodes, and attempt to explore it
		choice = rand() % choices.size();
		switch(choices[choice])
		{
		case 1: //Move Up
			choices.erase(choices.begin() + choice);  //Remove this from the choices list so that we don't attempt to explore it again
			if (curPos.y < MAZE_Y_SIZE - 1 && !maze[curPos.y + 1][curPos.x].isVisited) //If within maze bounds and if not an explored node
			{
				newPos.x = curPos.x;
				newPos.y = curPos.y + 1;
				maze[curPos.y][curPos.x].topWall = false;
				maze[newPos.y][newPos.x].bottomWall = false;
				generateMaze(newPos);
			}
			break;
		case 2: //Move right
			choices.erase(choices.begin() + choice);
			if (curPos.x < MAZE_X_SIZE - 1 && !maze[curPos.y][curPos.x + 1].isVisited)
			{
				newPos.x = curPos.x + 1;
				newPos.y = curPos.y;
				maze[curPos.y][curPos.x].rightWall = false;
				maze[newPos.y][newPos.x].leftWall = false;
				generateMaze(newPos);
				
			}
			break;
		case 3: //Move down
			choices.erase(choices.begin() + choice);
			if (curPos.y > 0 && !maze[curPos.y - 1][curPos.x].isVisited)
			{
				newPos.x = curPos.x;
				newPos.y = curPos.y - 1;
				maze[curPos.y][curPos.x].bottomWall = false;
				maze[newPos.y][newPos.x].topWall = false;
				generateMaze(newPos);
			}
			
			break;
		case 4: //Move left
			choices.erase(choices.begin() + choice);
			if (curPos.x > 0 && !maze[curPos.y][curPos.x - 1].isVisited)
			{
				newPos.x = curPos.x - 1;
				newPos.y = curPos.y;
				maze[curPos.y][curPos.x].leftWall = false;
				maze[newPos.y][newPos.x].rightWall = false;
				generateMaze(newPos);
			}
			break;
		}
		
	}




}

//This function converts a maze array index-X-coordinate to an OpenGL coordinate in [1.0, 1.0]
double convertMazeXToGLCoord(double mazeX)
{
	return ((mazeX) /(MAZE_X_SIZE)) * 2 - 1;
}

//This function converts a maze array index-Ycoordinate to an OpenGL coordinate in [1.0, 1.0]
double convertMazeYToGLCoord(double mazeY)
{
	return ((mazeY) / (MAZE_Y_SIZE)) * 2 - 1;
}

//This function calls all necessary logic to generate a new random maze and reset the user's position and moves
void resetMaze()
{
	userPos = startPos;
	initMaze();
	Point startPoint;
	startPoint.x = startPos.x;
	startPoint.y = startPos.y;
	generateMaze(startPoint);
}

//This function moves the legs by a bit
//Parameter factor is applied to the amount of movement (1 or -1 is usually used)
void moveLegs(double factor)
{
	if (legForward)
	{
		legAngle += 1.0 * factor;
		
	}
	else
	{
		legAngle -= 1.0 * factor;
		
	}

	if (legAngle >= 30)
	{
		legAngle = 30;
		legForward = !legForward;
	}

	if (legAngle <= -30)
	{
		legAngle = -30;
		legForward = !legForward;
	}
	
}

//This function moves the arms by a bit
//Parameter factor is applied to the amount of movement
void moveArms(double factor)
{
	if (arm1Forward)
	{
		arm1Angle += 1.0 * factor;
		
	}
	else
	{
		arm1Angle -= 1.0 * factor;
		
	}

	if (arm1Angle >= 45)
	{
		arm1Angle = 45;
		arm1Forward = !arm1Forward;
	}

	if (arm1Angle <= -45)
	{
		arm1Angle = -45;
		arm1Forward = !arm1Forward;
	}

	if (arm2Forward)
	{
		arm2Angle += 1.0 * factor;
	}
	else
	{
		arm2Angle -= 1.0 * factor;
	}

	if (arm2Angle >= 25)
	{
		arm2Angle = 25;
		arm2Forward = !arm2Forward;
	}

	if (arm2Angle <= -25)
	{
		arm2Angle = -25;
		arm2Forward = !arm2Forward;
	}
}

//This function is responsible for processing key presses through GLUT
//inputs -
//key - the key value to be processed
//val1 and val2 - not used
void keyBoardHandler (unsigned char key, int val1, int val2)
{
	vec4 displacement(0.0f, -0.1f, 0.0f, 1.0f);  //Unorientated the robot faces south

	mat4 rotationMatrix = rotate(mat4(1.0f), (float)angle, vec3(0.0f, 0.0f, 1.0f)); 

	displacement = rotationMatrix * displacement;

	switch(key)
	{
	case 'r': //Turn right 90 degrees
	case 'R':
		angle -= 90;
		glutPostRedisplay();
		break;
	case 'l': //Turn left 90 degrees
	case 'L':
		angle += 90;
		glutPostRedisplay();
		break;
	case 'f': //Move forwards
	case 'F':
		if (userPos.x + displacement[0] <= MAZE_X_SIZE && userPos.y + displacement[1] <= MAZE_Y_SIZE && userPos.x + displacement[0] >= 0 && userPos.y + displacement[1] >= 0)
		{
			userPos.x += displacement[0];
			userPos.y += displacement[1];
			moveLegs(1);
			moveArms(1);
			glutPostRedisplay();
		}
		break;
	case 'b': //Move backwards
	case 'B':
		if (userPos.x + displacement[0] <= MAZE_X_SIZE && userPos.y + displacement[1] <= MAZE_Y_SIZE && userPos.x + displacement[0] >= 0 && userPos.y + displacement[1] >= 0)
		{
			userPos.x -= displacement[0];
			userPos.y -= displacement[1];
			moveLegs(-1);
			moveArms(-1);
			glutPostRedisplay();
		}
		break;
	case 'm': //Reset Maze - new random maze
	case 'M':
		resetMaze();
		glutPostRedisplay();
		break;
	case 'c': //Return user to start of maze, but do not regenerate maze structure
	case 'C':
		userPos = startPos;
		glutPostRedisplay();
		break;
	case 'W': //Toggle wireframe mode
	case 'w':
		fillMode = !fillMode;
		if (!fillMode)
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
		glutPostRedisplay();
		break;
	case 'A': //Toggle high ambience mode
	case 'a':
		ambientMode = !ambientMode;
		glutPostRedisplay();
		break;
	}

}



//This function draws the robot
//Parameter totalMatrix is the transformation to be applied
void drawRobot(mat4 & projMatrix, mat4 & modelViewMatrix)
{
	mat4 robotMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(userPos.x), height/2, -convertMazeYToGLCoord(userPos.y)));
	robotMatrix = scale(robotMatrix, vec3(0.8f, 0.8f, 0.8f)); //80% of size of one maze cell
	robotMatrix = scale(robotMatrix, vec3(0.4f, 0.4f, 0.4f));  //It was easier to draw the robot higher than needed then scale him down

	//Robot overall orientation
	robotMatrix = rotate(robotMatrix, (float)angle, vec3(0.0f, 1.0f, 0.0f));

	//Left leg
	mat4 leftLegMatrix = robotMatrix;
	leftLegMatrix = translate(leftLegMatrix, vec3(-width/3, 0.0f, 0.0f));
	leftLegMatrix = translate(leftLegMatrix, vec3(0.0f, 0, 0.0f));
	leftLegMatrix = rotate(leftLegMatrix, -180.0f, vec3(1.0, 0.0f, 0.0f));
	leftLegMatrix = rotate(leftLegMatrix, (float)legAngle, vec3(1.0f, 0.0f, 0.0f)); //Left leg rotation angle
	leftLegMatrix = scale(leftLegMatrix, vec3(0.2f, 1.0f, 0.2f));

	greenCube.draw(projMatrix, leftLegMatrix);

	//Right leg
	mat4 rightLegMatrix = robotMatrix;
	rightLegMatrix = translate(rightLegMatrix, vec3(width/3, 0.0f, 0.0f));
	rightLegMatrix = translate(rightLegMatrix, vec3(0.0f, 0, 0.0f));
	rightLegMatrix = rotate(rightLegMatrix, -180.0f, vec3(1.0, 0.0f, 0.0f));
	rightLegMatrix = rotate(rightLegMatrix, (float)-legAngle, vec3(1.0f, 0.0f, 0.0f)); //Right leg rotation angle
	rightLegMatrix = scale(rightLegMatrix, vec3(0.2f, 1.0f, 0.2f));

	greenCube.draw(projMatrix, rightLegMatrix);

	mat4 bodyMatrix = translate(robotMatrix, vec3(0.0f, height, 0.0f));
	bodyMatrix = scale(bodyMatrix, vec3(0.9f, 1.0f, 0.9f));
	greenCube.draw(projMatrix, bodyMatrix);

	//Right arm top
	mat4 armMatrix1 = translate(robotMatrix, vec3(-width/4, height * 1.65, 0.0f));
	armMatrix1 = rotate(armMatrix1, 180.0f, vec3(0.0f, 1.0f, 0.0f));
	armMatrix1 = rotate(armMatrix1, (float)arm1Angle, vec3(0.0f, 1.0f, 0.0f)); //rotation angle
	armMatrix1 = translate(armMatrix1, vec3(width/2, 0.0f, 0.0f));
	matrixStack.push(armMatrix1); //Do not use  this scaling at next level of hiearchy
	armMatrix1 = scale(armMatrix1, vec3(0.7f, 0.2, 0.2f));
	greenCube.draw(projMatrix, armMatrix1);

	//Right arm bottom
	mat4 armMatrix2 = matrixStack.top();
	matrixStack.pop();
	armMatrix2 = translate(armMatrix2, vec3(width * 0.17, height/16, 0.0f));
	armMatrix2 = rotate(armMatrix2, 270.0f, vec3(0.0f, 0.0f, 1.0f));
	armMatrix2 = rotate(armMatrix2, (float)arm2Angle, vec3(0.0f, 0.0f, 1.0f));
	armMatrix2 = translate(armMatrix2, vec3(width/4, 0.0f, 0.0f));
	
	armMatrix2 = scale(armMatrix2, vec3(0.55f, 0.15, 0.10f));
	greenCube.draw(projMatrix, armMatrix2);

	/////////////////////////////////////////

	//Left arm top new
	armMatrix1 = translate(robotMatrix, vec3(width/4, height * 1.65, 0.0f));
	armMatrix1 = rotate(armMatrix1, 180.0f, vec3(0.0f, 1.0f, 0.0f));
	armMatrix1 = rotate(armMatrix1, (float)arm1Angle, vec3(0.0f, 1.0f, 0.0f)); //rotation angle
	armMatrix1 = translate(armMatrix1, vec3(-width/2, 0.0f, 0.0f));
	matrixStack.push(armMatrix1); //Do not use  this scaling at next level of hiearchy
	armMatrix1 = scale(armMatrix1, vec3(0.7f, 0.2, 0.2f));
	greenCube.draw(projMatrix, armMatrix1);

	//Left arm bottom new
	armMatrix2 = matrixStack.top();
	matrixStack.pop();
	armMatrix2 = translate(armMatrix2, vec3(-width * 0.17, height/16, 0.0f));
	armMatrix2 = rotate(armMatrix2, -270.0f, vec3(0.0f, 0.0f, 1.0f));
	armMatrix2 = rotate(armMatrix2, (float)arm2Angle, vec3(0.0f, 0.0f, 1.0f));
	armMatrix2 = translate(armMatrix2, vec3(-width/4, 0.0f, 0.0f));
	
	armMatrix2 = scale(armMatrix2, vec3(0.55f, 0.15, 0.10f));
	greenCube.draw(projMatrix, armMatrix2);

	//////////////////////////////////////////////

	//head
	mat4 headMatrix = translate(robotMatrix, vec3(0.0f, height*2.2, height/16));
	headMatrix = scale(headMatrix, vec3(0.6f, 0.6f, 0.8f));
	greenCube.draw(projMatrix, headMatrix);


}





//This function draws a white-green sphere on a stand (floating in the air)
//x/y/z - location of sphere-stand
//radius - radius of sphere component
//increment - increment used to find quads for sphere
//totalMatrix - transformation to be applied
void drawSphereOnStand(mat4 & projMatrix, mat4 & modelViewMatrix)
{

	mat4 newMatrix = translate(modelViewMatrix, vec3(0.0f, height*1.5, 0.0f));

	mat4 cubeMatrix = scale(newMatrix, vec3(0.75f, 0.5f, 0.75f));
	greenCube2.draw(projMatrix, cubeMatrix);

	newMatrix = translate(newMatrix, vec3(0.0f, height/6, 0.0f));
	standSphere.draw(projMatrix, newMatrix);
}


//Render function (AKA display function)
void render2()
{
	/////////////////
	//Deformations code
	
	//Update Logic
	double timeElapsed;

	switch (whichMethod)
	{
	case 1:					//Stanford Method
		timeElapsed = 0.00225;
		break;
	case 2:					//Georgia Institute Method
		timeElapsed = 0.00225;
		break;
	case 3:					//Non Linear Paper Method
		timeElapsed = 0.00225;
		break;
	default:
		timeElapsed = 0.005;
		break;
	}

	for (int i = 0; i < 10; i++)
	{				
		particleSystem -> doUpdate(timeElapsed);
	}

	//startTime = glutGet(GLUT_ELAPSED_TIME);
		
	particleSystem -> calculateNormals();
	///////////////////end of deformations code///////////////

	glClearColor(0.0, 0,0, 0.0);
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 
	mat4 projMatrix = perspective(45.0f, (float)aspectRatio, 0.001f, 100.0f); //Projection Matrix

	double cameraHeight = 0.5;  //Max height of camera if looking straight down at character

	//Rotate/scale the camera's eye point based on the mouse rotation/scale settings
	vec4 eyeToInterestVector(0, cameraHeight, 0, 1); //Displacement of eye position versus interest point
	
	mat4 transformMatrix = rotate(mat4(1.0f), x_angle, vec3(0.0f, 1.0f, 0.0f));
	transformMatrix = rotate(transformMatrix, y_angle, vec3(1.0f, 0.0f, 0.0f));
	transformMatrix = scale(transformMatrix, vec3(1.0f, scale_size, 1.0f));

	eyeToInterestVector = transformMatrix * eyeToInterestVector;

	vec3 eyeToInterestVector3D(eyeToInterestVector);

	//Have the camera look directly at the robot character
	//Rotate around him, or zoom in / out - based on the mouse camera settings
	vec3 eyePos(convertMazeXToGLCoord(userPos.x) + eyeToInterestVector3D[0], eyeToInterestVector3D[1], -convertMazeYToGLCoord(userPos.y) + eyeToInterestVector3D[2]);
	vec3 interestPoint(convertMazeXToGLCoord(userPos.x), 0.0, -(convertMazeYToGLCoord(userPos.y)));
	mat4 viewingMatrix = lookAt(eyePos, interestPoint, vec3(0.0, 1.0, 0.0));

	mat4 modelingMatrix = mat4(1.0f); //Identity matrix
	mat4 modelViewMatrix = viewingMatrix * modelingMatrix;


	
	//Use the pythagoeran theorem to calculate a radius for the sphere corresponding to the cell dimensions
	double radius = std::sqrt(width*width + height*height);

	//Goal cylinder
	mat4 shiftedMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(endPos.x) + width/2, height/1.25, -convertMazeYToGLCoord(endPos.y) - height/2));
	shiftedMatrix = scale(shiftedMatrix, vec3(0.7, 1.0, 0.7));
	cylinder.draw(projMatrix, shiftedMatrix);

	

	//A decoration at each corner
	shiftedMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(0) + width/2, 0, -convertMazeYToGLCoord(0) - height/2));
	drawSphereOnStand(projMatrix, shiftedMatrix);

	shiftedMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(MAZE_X_SIZE - 1) + width/2, 0, -convertMazeYToGLCoord(0) - height/2));
	drawSphereOnStand(projMatrix, shiftedMatrix);

	
	shiftedMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(0) + width/2, 0, -convertMazeYToGLCoord(MAZE_Y_SIZE - 1) - height/2));
	drawSphereOnStand(projMatrix, shiftedMatrix);
	
	
	shiftedMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(MAZE_X_SIZE - 1) + width/2, 0, -convertMazeYToGLCoord(MAZE_Y_SIZE - 1) - height/2));
	drawSphereOnStand(projMatrix, shiftedMatrix);
	
	float scaleFactor1 = 0.04; //Factor by which to shrink a cube along the needed dimension to make it into a wall
	float scaleFactor2 = 0.96; //Factor by which to shrink a cube along the other dimension so as to prevent z fighting amongst walls

	//Draw Maze
	for (int i = 0; i < MAZE_Y_SIZE; i++) //z dimension in 3d
	{
		for (int j = 0; j < MAZE_X_SIZE; j++) //x dimension in 3d
		{

			if (maze[i][j].topWall)
			{
				
				mat4 wallMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(j) + width / 2, yHeight, -convertMazeYToGLCoord(i+1)));
				wallMatrix = scale(wallMatrix, vec3(1.0f, 1.0f, scaleFactor1));
				blueCube.draw(projMatrix, wallMatrix);
				
				

				
			}
			if (maze[i][j].rightWall)
			{
				

				mat4 wallMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(j+1), yHeight, -convertMazeYToGLCoord(i) - height/2));
				wallMatrix = scale(wallMatrix, vec3(scaleFactor1, 1.0f, scaleFactor2));  //1st parameter should be 0.1
				blueCube.draw(projMatrix, wallMatrix);
				
				
				
			}

			if (maze[i][j].bottomWall)
			{

				mat4 wallMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(j) + width / 2, yHeight, -convertMazeYToGLCoord(i)));
				wallMatrix = scale(wallMatrix, vec3(1.0f, 1.0f, scaleFactor1));
				blueCube.draw(projMatrix, wallMatrix);
				
			}
			
			if (maze[i][j].leftWall)
			{

				mat4 wallMatrix = translate(modelViewMatrix, vec3(convertMazeXToGLCoord(j), yHeight, -convertMazeYToGLCoord(i) - height/2));
				wallMatrix = scale(wallMatrix, vec3(scaleFactor1, 1.0f, scaleFactor2));  //1st parameter should be 0.1
				blueCube.draw(projMatrix, wallMatrix);
				
				
			}
			
		}
	}
	

	modelViewMatrix = viewingMatrix * modelingMatrix;
	vec4 robotPos = vec4(convertMazeXToGLCoord(userPos.x), 0, -convertMazeYToGLCoord(userPos.y), 1);
	drawRobot(projMatrix, modelViewMatrix);
	
	////////////////////////////////////////
	//Deformations code
	const float SCALEFACTOR = 0.05f;
	mat4 tetraMatrix = scale(modelViewMatrix, vec3(SCALEFACTOR, SCALEFACTOR, SCALEFACTOR));


	particleSystem -> doRender(timeElapsed * 4, projMatrix, tetraMatrix);
	///////////////////////////////////////////

	glutSwapBuffers();
	
}

//Window resize GLUT event function
//Parameters windowWidth/windowHeight - new dimensions of the window
void resize2(int windowWidth, int windowHeight)
{
	aspectRatio = (float) windowWidth / windowHeight;
	glViewport(0.0f, 0.0f, windowWidth, windowHeight);
	glutPostRedisplay();
}


//Main function
int main()
{
	////////////////////
	//Deformations Code
	logger = new Logger();
	int vertexCount = 0;
	int tetraCount = 0;
	Vertex * vertexList = NULL;
	int * tetraList = NULL;
	TetraMeshReader theReader;
	//if (theReader.openFile("house2.node", "house2.ele"))
	if (theReader.openFile("chrisSimpler.node", "chrisSimpler.ele"))
	{
		bool loadSucceeded = theReader.loadData(vertexList, vertexCount, tetraList, tetraCount, logger);
		
		theReader.closeFile();

		if (loadSucceeded && vertexList != NULL && tetraList != NULL)
		{
			switch(whichMethod)
			{
			case 1:
				particleSystem = new StanfordSystem(vertexList, vertexCount, tetraList, tetraCount, logger);
				break;
			case 2:
				particleSystem = new GeorgiaInstituteSystem(vertexList, vertexCount, tetraList, tetraCount, logger);
				break;
			case 3:
				particleSystem = new NonlinearMethodSystem(vertexList, vertexCount, tetraList, tetraCount, logger);
				break;
			default:
				cerr << "Incorrect system identifier -- defaulting to stanford system" << endl;
				particleSystem = new StanfordSystem(vertexList, vertexCount, tetraList, tetraCount, logger);
				break;
			}

			//particleSystem -> loadSpecialState();
		}
		else
		{
			cerr << "Program cannot run with error in loading input data contents" << endl;
			system("pause");
			return 1;
		}
	}
	else
	{
		cerr << "Unable to execute program without input data" << endl;
		system("pause");
		return 1;
	}

	///////////////////

	srand(time(0));  //Seed random numbers based on current time
	
	//Set maze start and goal points
	startPos.x = 0.5;
	startPos.y = 0.5;
	endPos.x = MAZE_X_SIZE - 1;
	endPos.y = MAZE_Y_SIZE - 1;

	width = 2.0 / MAZE_X_SIZE;
	height = 2.0 / MAZE_Y_SIZE;
	yHeight = width / 2;

	resetMaze();

	

	//Glut initialization
	int mode = GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH;
	
	glutInitDisplayMode(mode);
	glutInitWindowSize(1000,800);
	glutCreateWindow("Maze Generator 3D");
	glutDisplayFunc(render2);
	glutIdleFunc(render2);
	glutKeyboardFunc(keyBoardHandler);
	glutReshapeFunc(resize2);
	glutMouseFunc(mymouse); 
    glutMotionFunc(mymotion);


	glEnable(GL_DEPTH_TEST);

	GLenum errorCode = glewInit();

	programObject = SetupGLSL("maze");

	double radius = std::sqrt(width*width + height*height);
	standSphere.init(radius / 5, 12, 12, vec3(0.5f, 1.0f, 0.5f), vec3(0.25f, 0.5f, 0.25f));
	standSphere.initVBO();
	

	cylinder.init(radius/2, radius/3, height, 10, 10, vec3(1.0f, 0.0f, 0.0f), vec3(0.5f, 0.0f, 0.0f));
	cylinder.initVBO();
	

	blueCube.init(width, yHeight, height, vec3(0.0f, 0.0f, 1.0f), vec3(0.0f, 0.0f, 1.0f));
	blueCube.initVBO();

	greenCube.init(width * 1.0, yHeight * 3, height, vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f));
	greenCube.initVBO();
	greenCube.setShininess(10000);

	greenCube2.init(width, yHeight, height, vec3(0.0f, 1.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f));
	greenCube2.initVBO();
	greenCube2.setShininess(5000);

	//////////////Deformations code
	particleSystem->initVBOs();
	particleSystem->setProgramObject(programObject);
	/////////////////////////////////

	glutMainLoop();

	/////////////////
	//Deformations code
	delete particleSystem;
	delete logger;
	////////////////
	
	return 0;
}
