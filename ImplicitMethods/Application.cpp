//Cloth simulation using implicit methods to implement springs for a grid of particles
//Features:
//	*Implements springs between particles defined in an edge list using hooke's law with a damping component
//	*Has an earth gravity component to drag the cloth downwards
//	*Uses implicit methods, allowing for a higher spring constant (variable ks in class ParticleSystem)
//	*Uses a conjugate gradient solver to implement implicit methods
//	*Gourad shading for lighting
//Controls:
//Keyboard:
//	SPACE: Begin animating
//	A: Decrease rest length between particles
//	D: Increase rest length between particles
//	W: Increase earth gravity (starts at 9.8 m/sec)
//	S: Decrease earth gravity
//	Z: toggle wireframe mode
//	X: toggle informational text display
//  E: run an explicit implementation of the simulation (useful for comparison; most obvious if you turn automatic implicit animation off with space bar)
//  R: reset the simulation
//  I: render to a series of numbered images so that they can be combined into a video
//	O: toggle logging of positions/velocities of constrained particles (only if DEBUGGING macro is #defined in Logger.h)
//	P: toggle complete logging (only if DEBUGGING macro is #defined in Logger.h)
//Mouse:
//	Left button - hold this whie dragging the mouse to change the rotation angle of the piece of cloth shown
//	Middle button - zoom in
//	Right button - zoom out
//Written by Chris Jacobsen with advisement from Professor Huamin Wang

//Based on:
//Deformation Papers:
//http://www.math.ucla.edu/~jteran/papers/TSNF03.pdf – Finite Volume Methods for the Simulation of Skeletal Muscle
//By R. Fedkiw et. al
//Deformation Method #1 - Class Stanford System
//http://graphics.berkeley.edu/papers/Obrien-GMA-1999-08/Obrien-GMA-1999-08.pdf – Graphical Modeling and Animation of Brittle Fracture
//By James O'Brien and Jessica Hodgkins
//Deformation Method #2 - Class GeorgiaInstituteSystem
//http://www-ljk.imag.fr/Publications/Basilic/com.lmc.publi.PUBLI_Article@11f6a0378d9_18c74/tensile.pdf – Simple, yet Accurate Nonlinear Tensile Stiffness
//Pascal Volino et. al
//Deformation Method #3 - Class NonlinearMethodSystem

#include <gl/glew.h>

#include <gl/glut.h>
#include <gl/GLU.H>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>

#include <string>
#include <sstream>
#include <iostream>
#include "StanfordSystem.h"
#include "GeorgiaInstituteSystem.h"
#include "NonlinearMethodSystem.h"
#include "ViewManager.h"
#include "Keyboard.h"
#include "TetraMeshReader.h"

using namespace std;

GLuint SetupGLSL(char *fileName);

//Note: the reason these were declared globally is to accommodate Glut's function calling system
ParticleSystem * particleSystem;	//The main particle system
ViewManager viewManager;			//Instance of the view manager to allow user view control
Keyboard * keyboard;				//Instance of the Keyboard class to process key presses
Logger * logger;					//Instance of Logger class to perform all logging
const int whichMethod = 1;			//1 for stanford method.  2 for georgia Institute Method.  3 for NonLinear Paper method.
const int whichModel = 4;

double ar = 0;

int summedTime = 0;
int frameCount = 0;
const int FRAME_COUNT_LIMIT = 1000;

//This function is called for rendering by GLUT
void render()
{
	//Timing mechanism for performance evaluation
	int startTime, endTime;
	startTime = glutGet(GLUT_ELAPSED_TIME);
	
	//Update Logic
	double timeElapsed;

	if (whichModel == 3)
	{
		//Gigantic dragon model -- need bigger constants
		switch (whichMethod)
		{
		case 1:					//Stanford Method
			timeElapsed = 0.00225;
			particleSystem->setConstants(2800,2800);
			break;
		case 2:					//Georgia Institute Method
			timeElapsed = 0.00225;
			particleSystem->setConstants(2800,2800);
			break;
		case 3:					//Non Linear Paper Method
			timeElapsed = 0.00225;
			particleSystem->setConstants(7000,7000);
			break;
		default:
			timeElapsed = 0.005;
			break;
		}
	}
	else if (whichModel == 1 || whichModel == 2)
	{
		//Regular model of some sort - do not need quite as big of constants
		switch (whichMethod)
		{
		case 1:					//Stanford Method
			timeElapsed = 0.00225;
			particleSystem->setConstants(700,700);
			break;
		case 2:					//Georgia Institute Method
			timeElapsed = 0.00225;
			particleSystem->setConstants(700,700);
			break;
		case 3:					//Non Linear Paper Method
			timeElapsed = 0.00225;
			particleSystem->setConstants(600,600);
			break;
		default:
			timeElapsed = 0.005;
			break;
		}
	}
	else 
	{
		//The 'simpler' tetrahedral model
		switch (whichMethod)
		{
		case 1:					//Stanford Method
			timeElapsed = 0.005;
			particleSystem->setConstants(700,700);
			break;
		case 2:					//Georgia Institute Method
			timeElapsed = 0.00225;
			particleSystem->setConstants(700,700);
			break;
		case 3:					//Non Linear Paper Method
			timeElapsed = 0.00225;
			particleSystem->setConstants(600,600);
			break;
		default:
			timeElapsed = 0.005;
			break;
		}
	}
	
	//for (int i = 0; i < 1; i++)
	for (int i = 0; i < 10; i++)
	{				
		particleSystem -> doUpdate(timeElapsed);
	}
	
	//startTime = glutGet(GLUT_ELAPSED_TIME);
		
	particleSystem -> calculateNormals();

	viewManager.doUpdate(timeElapsed);

	//Render Logic
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glm::mat4 projMatrix = glm::perspective(45.0f, (float)ar, 0.1f, 100.0f); //Projection Matrix
	double cameraHeight = 0.5;  //Max height of camera if looking straight down at character
	
	glm::mat4 floorModelViewMatrix = viewManager.doTransform();
	glm::mat4 tetraModelViewMatrix = floorModelViewMatrix;

	if (whichModel == 2)
	{
		tetraModelViewMatrix= glm::translate(tetraModelViewMatrix, glm::vec3(-10.0f, 0, 7.0f));
	}

	if (whichModel == 3)
	{
		tetraModelViewMatrix= glm::scale(tetraModelViewMatrix, glm::vec3(0.25, 0.25, 0.25));
		tetraModelViewMatrix= glm::translate(tetraModelViewMatrix, glm::vec3(-45.0f, -12, 30.0f));

	}

	particleSystem -> doRender(timeElapsed * 4, projMatrix, floorModelViewMatrix, tetraModelViewMatrix);

	glFlush();

	glutSwapBuffers();

	endTime = glutGet(GLUT_ELAPSED_TIME);
	//cout << "Total time for frame was: " << (double)(endTime - startTime) / 1000 << " seconds" << endl;
	frameCount++;
	summedTime += (endTime - startTime);

	if (frameCount == FRAME_COUNT_LIMIT) 
	{
		double fps = frameCount / (summedTime / double(1000));
		cout << "Average FPS for " << FRAME_COUNT_LIMIT << " frames is: " << fps << endl;
		frameCount = 0;
		summedTime = 0;
 	}

}

//This function processes window resize events, ensuring the aspect ratio and anything else dependent on window dimensions is correct
//Parameters width and height are the dimensions of the window
void resize(int width, int height)
{
	//Minimizing the window results in sizes of 0, which causes problems with video generation.
	//Do not allow the width or height to be set to 0.
	if (width == 0 || height == 0)
	{
		return;
	}

	//Calculate the aspect ratio using the x and y dimensions
	ar = (double) width / height;

	//Set the Open GL viewport to the window's dimensions
	glViewport(0, 0, width, height);

	//Go back to MODEL VIEW matrix mode
	particleSystem -> setWindowDimensions(width,height);
}

//This function processes mouse clicks (when the button is pressed)
//Parameters:
//button - which button was pressed
//state - was it up or down?
//x - x coordinate of click
//y - y coordinte of click
void mouseClick(int button, int state, int x, int y)
{
	viewManager.mouseClick(button,state,x,y);
}

//This function processes mouse motion.
//It only takes action if a button is being held
//Parameters:
//x - x coordinate of movement
//y - y coordinate of movement
void mouseMove(int x, int y)
{
	viewManager.mouseMove(x,y);
}

//This function processes keyPress events (a keyboard button going down)
//Parameter key - the button being pressed
//The other parameters are not used
void keyPressed (unsigned char key, int mystery, int mystery2)
{
	keyboard -> keyPressed(key);
}

//This function processes keyRelease events (a keyboard button going up)
//Parameter key - the button being released
//The other parameters are not used
void keyReleased (unsigned char key, int mystery, int mystery2)
{
	keyboard -> keyReleased(key);
}


//Main function
int main(int argCount, char **argValue)
{
	logger = new Logger();

	int vertexCount = 0;
	int tetraCount = 0;
	Vertex * vertexList = NULL;
	int * tetraList = NULL;
	TetraMeshReader theReader;
	
	bool loadSucceeded = false;
	if (whichModel == 1)
	{
		loadSucceeded = theReader.openFile("house2.node", "house2.ele");
	}
	else if (whichModel == 2)
	{
		loadSucceeded = theReader.openFile("P.node", "P.ele");
	}
	else if (whichModel == 3)
	{
		loadSucceeded = theReader.openFile("dragon.node", "dragon.ele");
	}
	else
	{
		loadSucceeded = theReader.openFile("chrisSimpler.node", "chrisSimpler.ele");
	}

	if (loadSucceeded)
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
			
			keyboard = new Keyboard(particleSystem, &viewManager, logger);

			glutInit(&argCount,argValue);
			glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE | GLUT_RGBA);

			glutInitWindowPosition(0,0);
			int windowWidth = 1000;
			int windowHeight = 700;
			glutInitWindowSize(windowWidth, windowHeight);
			resize(windowWidth, windowHeight);
			particleSystem -> setWindowDimensions(windowWidth, windowHeight);
			glutCreateWindow("Implicit Methods");

			//Note: this must be AFTER the create window function call to work properly
			glEnable(GL_DEPTH_TEST);
			glDepthFunc(GL_LEQUAL);
			glShadeModel(GL_SMOOTH);

			GLenum errorCode = glewInit();
			GLuint programObject = SetupGLSL("maze");
			particleSystem->initVBOs();
			particleSystem->setProgramObject(programObject);

			glClearColor(0.0f,0.0f,0.0f,0.0f);
	
			glutDisplayFunc(render);
			glutIdleFunc(render);
			glutReshapeFunc(resize);
			glutMouseFunc(mouseClick);
			glutMotionFunc(mouseMove);
			glutKeyboardFunc(keyPressed);
			glutKeyboardUpFunc(keyReleased);
			
			glutMainLoop();

			delete keyboard;
			delete particleSystem;
		}
		else
		{
			cerr << "Program cannot run with error in loading input data contents" << endl;
			system("pause");
		}
	}
	else
	{
		cerr << "Unable to execute program without input data" << endl;
		system("pause");
	}

	delete logger;

	return 0;

}
