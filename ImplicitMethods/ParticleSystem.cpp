#include <gl/glew.h>

#include <gl/glut.h>
#include <gl/GLU.h>
#include <cmath>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include "targa.h"
#include "Macros.h"

#include "ParticleSystem.h"



using namespace std;

////////////////////////////////////////
//Main ParticleSystem class implementation
////It builds a grid of particles (represented by Particle objects)
//These are then held together by springs using hooke's laws.  One edge corresponds to each spring.
//Hooke's law is used to represent the springs, with both a spring component and a damping comopnent
//Implicit methods are used for the integration.  This requires a solving system but allows much larger spring constants without instability
//and potentially allows more efficient implementation
//A gravity component is also present
////////////////////////////////////////
//Issue conclusions - for 3 tetrahedra, it has issues in C++.  It does not seem to have them in Matlab.
//It seems to "blast off" upwards.  It does this without collisions response even present.
//For a couple iterations, the computed def vertices were the same between C++ and Matlab
//I think we need to run the implementations side by side, only printing defVertices, and finding out when and where they divege
//Then, at that iteration, break it doesn by strain, force, etc, and find out the big piece that's wrong
//Then drill down further to figure out why.

//Based on:
//http://graphics.snu.ac.kr/~kjchoi/publication/cloth.pdf - Explicit and implicit formulas for hooke's law - page 3
//http://www.amath.unc.edu/Faculty/mucha/Reprints/SCAclothcontrolpreprint.pdf - Formulas for A and b in Ax = b system - page 5
//http://en.wikipedia.org/wiki/Conjugate_gradient - Algorithm for conjugate gradient (solves Ax = b)

const double epsilon = 1e-12;	//Used to check approximate equality to 0
extern const int DIMENSION;		//DIMENSION of system (3 for 3D)

//Constructor - initializes particles and settings
ParticleSystem::ParticleSystem(Vertex * vertexList, int vertexCount, int * tetraList, int tetraCount, Logger * logger)
{
	this -> logger = logger;
	dimensionSquared = DIMENSION * DIMENSION;

	/////////////////////////////////////
	double K = 55 * 2;					//Bulk Modulus
	mu = 55 * 2;						//Shear modulus (Lame's second parameter)
	lambda = K - (2.0/3) * mu;			//Lame's first parameter
	//kd = 1.93;
	kd = 7;

	restitution = 0.5;	//good
	//restitution = 1;	//jitters permanently
	ud = 0.3;	//good
	//ud = 0.9;	//shoves to left
	us = 0;

	showNormals = false;
	showOrg = 1;
	
	//numVertices = 4;
	numVertices = vertexCount;
	
	//orgVertices = new Vertex[numVertices];
	defVertices = new Vertex[numVertices];
	orgVertices = vertexList;
	
	const double height = 1.0;
	//const double height = -3.0;
	const double triangleDepth = sqrt(6.75);
	
	//orgVertices[0].position[0] = -1.5;         orgVertices[1].position[0] =   0;                        orgVertices[2].position[0] = 1.5;           orgVertices[3].position[0] =  0;
    //orgVertices[0].position[1] = 0 + height;   orgVertices[1].position[1] =   triangleDepth + height;   orgVertices[2].position[1] = 0 + height;    orgVertices[3].position[1] = 0 + height; 
    //orgVertices[0].position[2] =  0;           orgVertices[1].position[2] =   -triangleDepth / 2;       orgVertices[2].position[2] =  0;            orgVertices[3].position[2] = -triangleDepth;

	reset();

	//numTetra = 1;
	//tetraList = new int[numTetra * 4];
	numTetra = tetraCount;
	this -> tetraList = tetraList;
	
	//Note: these are in CLOCKWISE ORDER.  Winding must be consistent.
	//tetraList[0] = 0;
	//tetraList[1] = 1;
	//tetraList[2] = 2;
	//tetraList[3] = 3;
   		
	earthGravityValue = 9.8;	//meters per second for earth gravity force
	//Set number of rows and columns for grid here
		
	//Compute half of width/height - allows grid to be relatively centered while rendering
	//halfWidth = sheetWidth / 2;
	//halfHeight = sheetHeight / 2;

	halfWidth = 10; //temp!
	halfHeight = 10; //temp!
		
	//reset();  //Set up the particle positions / velocities

	normals = new double[DIMENSION * numTetra * 4];

	//vertexNormals = new double [(DIMENSION + 1) * numVertices];
	tetraCounts = new int [numVertices];

	massMatrix = new double [numVertices];
	currentForce = new double[numVertices * DIMENSION];
	zeroVector = new double [numVertices * DIMENSION];

	//Allocate buffers to hold portions of the screen grabbed / swapped for video generation
	//We allocate enough space to work with the entire screen because the window size must be less than the screen
	int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
	int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
	screenBuffer = new uint8_t[screenWidth * screenHeight * 4];
	screenRowTemp = new uint8_t [screenWidth * 4];
	
	
	
	//Create mass matrix.  Represent it with a vector, since it's a diagonal matrix.
	//Referred to http://en.wikipedia.org/wiki/Mass_matrix for overview of mass matrices
	for (int i = 0; i < numVertices; i++)
	{
		massMatrix[i] = 1;				
	}
	
	for (int i = 0; i < DIMENSION * numVertices; i++)
	{
		zeroVector[i] = 0;
	}
	
	//Set values for boolean variables
	showInfoText = true;
	
	isAnimating = false;
	renderMode = 1;
	renderToImage = false;
	frameNumber = 1;
	timeSinceVideoWrite = 0.0;

	

	for (int i = 0; i < TEXT_SIZE; i++)
	{
		text[i] = '\0';
	}

	windowWidth = windowHeight = 0;

	system("del images\\*.tga");  //DOS command to remove all tga image files from previous executions of this program

	#ifdef DEBUGGING
	//logger -> loggingLevel = logger -> FULL;
	if (logger -> isLogging)
	{
		cout << "Version 5 - 1st explicit implementation of georgia institute method (brittle fracture paper)" << endl;
		//logger -> printParticles(particles, numParticles, logger -> MEDIUM);
		//logger -> printEdges(edgeList, numEdges, logger -> MEDIUM);
	}
	#endif

	////////////////////////////LIGHTING////////////////////////////////////////////////
	//If lighting is calculated in eye space, the eye position basically is the origin - use this for the default
	eyePos[0] = 0.0f;
	eyePos[1] = 0.0f;
	eyePos[2] = -5.0f;

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

	//Turning off specular for now unless it's needed (and verified as good looking)
	//lightSpecular[0] = 1;
	//lightSpecular[1] = 1;
	//lightSpecular[2] = 1;
	lightSpecular[0] = 0;
	lightSpecular[1] = 0;
	lightSpecular[2] = 0;
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

	lightColor[0] = 1.0f;
	lightColor[1] = 0.1f;
	lightColor[2] = 0.1f;
	lightColor[3] = 1.0f;

	matShininess[0] = 10000;
	
	ambientMode = false;
}

//Destructor - free all memory for dynamically allocated arrays
ParticleSystem::~ParticleSystem()
{
	

	delete [] orgVertices;
	delete [] defVertices;
	delete [] constraintParticles;
	
	
	delete [] massMatrix;
	delete [] zeroVector;
	delete [] currentForce;
	delete [] screenBuffer;
	delete [] screenRowTemp;

}

//This method resets the simulation by returning all particles to their original positions and velocities
//It is also used by the constructor to initialize all particles' positions and velocities, forming a sheet (grid) of cloth
//The particles essentially correspond to vertices in a graph
//This method does NOT change view control, changes to the gravity setting, or changes to the rest length setting
void ParticleSystem::reset()
{
	//for (int i = 0; i < rows; i++)
	//{
	//	for (int j = 0; j < cols; j++)
	//	{
	//		particles[i * cols + j].position[0] = 0 + j * restLengthValues[0];
	//		particles[i * cols + j].position[2] = 0 + -i * restLengthValues[1];
	//		particles[i * cols + j].position[1] = 0;
	//		//particles[i * cols + j].position.y = 0 + -i * padding;		//Useful if a vertical sheet is needed
	//		//particles[i * cols + j].position.z = 0;						//Useful if a vertical sheet is needed
	//		for (int k = 0; k < DIMENSION; k++)
	//		{
	//			particles[i * cols + j].velocity[k] = 0;
	//		}

	//	}
	//}

	for (int i = 0; i < numVertices; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			defVertices[i].position[j] = orgVertices[i].position[j];
			//orgVertices[i].velocity[j] = defVertices[i].velocity[j] = (j == 0?3:0);
			orgVertices[i].velocity[j] = defVertices[i].velocity[j] = (j == 0?0:0);
		}
	}

	doTransform();

}

//Update Method - Implements one time step for the animation
//Uses hooke's law with implicit methods to construct matrix A and vector b then calls the conjugate gradient method to find deltaV
//This is then used to update the particle velocities and in turn the positions
//Parameter - deltaT - Amount of time elapsed to use in integrating.  Type double. 
void ParticleSystem::doUpdate(double deltaT)
{
	
}

void ParticleSystem::doCollisionDetectionAndResponse(double deltaT)
{
	//for k = 1:size(defVertices,2)
	for (int i = 0; i < numVertices; i++)
		{
			//if defVertices(2,k) < floor
			///if (defVertices[i].position[1] < 1)  //TEMP HACK
			if (defVertices[i].position[1] < -4)
			{
				#ifdef DEBUGGING
				if (logger -> isLogging)
				{
					logger ->printVelocitiesAndPositions(defVertices,numVertices,"Positions during collision response", "Velocities", logger->MEDIUM);
				}
				#endif
				
				//%From: http://gafferongames.com/virtualgo/collision-response-and-coulomb-friction/
				//j = max(-(1 + restitution) * dot(inVelocities(:,k), [0 1 0]'), 0); %Magnitude of reaction force without rotational inertia

				double verticalNormal[3] = {0, 1, 0};

				double dotProduct = 0;
				for (int j = 0; j < DIMENSION; j++)
				{
					dotProduct += defVertices[i].velocity[j] * verticalNormal[j];
				}
				
				double jr = -(1 + restitution) * dotProduct;
				if (jr < 0)
				{
					jr = 0;
				}

				#ifdef DEBUGGING
				if (logger -> isLogging)
				{
					cout << "jr is: " << jr << endl;
				}
				#endif

				//js = us * j;
				//jd = ud * j;
				double js = us * jr;
				double jd = ud * jr;
                
				//tangent = inVelocities(:,k) - dot(inVelocities(:,k), [0 1 0]') * [0 1 0]';
				double tangent[3];
				double magnitude = 0; //This is the magnitude of the tangent vector
				for (int j = 0; j < DIMENSION; j++)
				{
					tangent[j] = defVertices[i].velocity[j] - dotProduct * verticalNormal[j];
					magnitude += tangent[j] * tangent[j];
				}
				magnitude = sqrt(magnitude);

				//if norm(tangent) ~= 0 && dot(inVelocities(:,k),[0 1 0]' ~= 0);
				//	tangent = tangent / norm(tangent);       
				//else %Note: decided to omit the "sum of external forces" case for now - just setting it to zero in this case
				//	tangent = [0 0 0]'; %Being really explicit in the code for now
				//end
				if (magnitude > epsilon && abs(dotProduct) > epsilon) //STABILITY FIX - dot product check
				{
					for (int j = 0; j < DIMENSION; j++)
					{
						tangent[j] /= magnitude;
					}
				}
				else //Precision problem - magnitude rounded to zero.  We cannot normalize the vector with it, or we'll get division by zero.
				{
					tangent[0] = 0;		//STABILITY FIX - 0 vector instead of 1 0 0 vector
					tangent[1] = 0;
					tangent[2] = 0;
				}

				#ifdef DEBUGGING
				if (logger -> isLogging)
				{
					cout << "Tangent vector is: " << endl;
					for (int j = 0; j < DIMENSION; j++)
					{
						cout << tangent[j] << " ";
					}
					cout << endl;
				}
				#endif
				     
				double dotProduct2 = 0; //This is the dot product of the velocity and the tangent vectors
				for (int j = 0; j < DIMENSION; j++)
				{
					dotProduct2 += defVertices[i].velocity[j] * tangent[j];
				}

				//if dot(inVelocities(:,k),tangent) == 0 && dot(inVelocities(:,k),tangent) <= js
				//	jf =  -(dot(inVelocities(:,k), tangent)) * tangent;
				//else
				//	jf = -jd * tangent;
				//end
				double jf[3];
				if (abs(dotProduct2) < epsilon && dotProduct2 <= js) //STABILITY FIX - changed || to &&
				{
					for (int j = 0; j < DIMENSION; j++)
					{
						jf[j] = -dotProduct2 * tangent[j];
					}
				}
				else
				{
					for (int j = 0; j < DIMENSION; j++)
					{
						jf[j] = -jd * tangent[j];
					}
				}

				#ifdef DEBUGGING
				if (logger -> isLogging)
				{
					cout << "jf:" << endl;
					for (int j = 0; j < DIMENSION; j++)
					{
						cout << jf[j] << " ";
					}
					cout << endl;
			    }
				#endif
        
				//inVelocities(:,k) = inVelocities(:,k) + j * [0 1 0]' + jf; %Apply instantaneous reaction impulse force
				//defVertices(:,k) = defVertices(:,k) + (j * [0 1 0]' + jf) * elapsedTime;
				for (int j = 0; j < DIMENSION; j++)
				{
					defVertices[i].velocity[j] += jr * verticalNormal[j] + jf[j];
					defVertices[i].position[j] += (jr * verticalNormal[j] + jf[j]) * deltaT;
				}

				#ifdef DEBUGGING
				if (logger -> isLogging)
				{
					cout << "Final positions/velocities for your viewing pleasure..." << endl;
					logger ->printVelocitiesAndPositions(defVertices,numVertices,"Positions after collision response", "Velocities", logger->LIGHT);
				}
				#endif
			}
		//end
		} 


}



////Method to calculate the normals used for lighting
////Note that since the cross product is only defined in 3 dimensions, this method only works properly for 3 dimensions
void ParticleSystem::calculateNormals()
{
//	//Cross product and this function only work if DIMENSION == 3
//
	double vectorDifferenceA[DIMENSION];  //1st vector formed for each cross product
	double vectorDifferenceB[DIMENSION];  //2nd vector formed for each cross product
	double crossProductResult[DIMENSION]; //Holds the result of a cross product from the macro

	//Initialize items to 0
	for (int i = 0; i < numVertices; i++)
	{
		for (int j = 0; j < DIMENSION; j++)
		{
			normals[numVertices * j + i] = 0;
			defVertices[i].vertexNormal[j] = 0;
		}
		tetraCounts[i] = 0;
	}

	for (int currentTetrad = 0; currentTetrad < numTetra; currentTetrad++)
	{
		for (int i = 0; i < DIMENSION; i++)
		{
			vectorDifferenceA[i] = defVertices[tetraList[3 * numTetra + currentTetrad]].position[i] - defVertices[tetraList[1 * numTetra + currentTetrad]].position[i];
			vectorDifferenceB[i] = defVertices[tetraList[1 * numTetra + currentTetrad]].position[i] - defVertices[tetraList[0 * numTetra + currentTetrad]].position[i];
		}

		//crossProductGeneral(crossProductResult, vectorDifferenceA, vectorDifferenceB);
		crossProductGeneral(crossProductResult, vectorDifferenceA, vectorDifferenceB);

		for (int i = 0; i < DIMENSION; i++)
		{
			defVertices[tetraList[3 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			defVertices[tetraList[1 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			defVertices[tetraList[0 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[3 * numTetra + currentTetrad]] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[1 * numTetra + currentTetrad]] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[0 * numTetra + currentTetrad]] += crossProductResult[i];
		}
		tetraCounts[tetraList[3 * numTetra + currentTetrad]]++;
		tetraCounts[tetraList[1 * numTetra + currentTetrad]]++;
		tetraCounts[tetraList[0 * numTetra + currentTetrad]]++;

		//glVertex3f(defVertices[tetraList[3 * numTetra + i]].position[0], defVertices[tetraList[3 * numTetra + i]].position[1], defVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3
		//glVertex3f(defVertices[tetraList[1 * numTetra + i]].position[0], defVertices[tetraList[1 * numTetra + i]].position[1], defVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
		//glVertex3f(defVertices[tetraList[0 * numTetra + i]].position[0], defVertices[tetraList[0 * numTetra + i]].position[1], defVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0
		
		for (int i = 0; i < DIMENSION; i++)
		{
			vectorDifferenceA[i] = defVertices[tetraList[2 * numTetra + currentTetrad]].position[i] - defVertices[tetraList[1 * numTetra + currentTetrad]].position[i];
			vectorDifferenceB[i] = defVertices[tetraList[1 * numTetra + currentTetrad]].position[i] - defVertices[tetraList[3 * numTetra + currentTetrad]].position[i];
		}

		//crossProductGeneral(crossProductResult, vectorDifferenceA, vectorDifferenceB);
		crossProductGeneral(crossProductResult, vectorDifferenceA, vectorDifferenceB);

		for (int i = 0; i < DIMENSION; i++)
		{
			defVertices[tetraList[2 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			defVertices[tetraList[1 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			defVertices[tetraList[3 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[2 * numTetra + currentTetrad]] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[1 * numTetra + currentTetrad]] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[3 * numTetra + currentTetrad]] += crossProductResult[i];
		}
		tetraCounts[tetraList[2 * numTetra + currentTetrad]]++;
		tetraCounts[tetraList[1 * numTetra + currentTetrad]]++;
		tetraCounts[tetraList[3 * numTetra + currentTetrad]]++;

		//glVertex3f(defVertices[tetraList[2 * numTetra + i]].position[0], defVertices[tetraList[2 * numTetra + i]].position[1], defVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2
		//glVertex3f(defVertices[tetraList[1 * numTetra + i]].position[0], defVertices[tetraList[1 * numTetra + i]].position[1], defVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
		//glVertex3f(defVertices[tetraList[3 * numTetra + i]].position[0], defVertices[tetraList[3 * numTetra + i]].position[1], defVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3

		for (int i = 0; i < DIMENSION; i++)
		{
			vectorDifferenceA[i] =  defVertices[tetraList[2 * numTetra + currentTetrad]].position[i] - defVertices[tetraList[3 * numTetra + currentTetrad]].position[i];
			vectorDifferenceB[i] =  defVertices[tetraList[3 * numTetra + currentTetrad]].position[i] - defVertices[tetraList[0 * numTetra + currentTetrad]].position[i];
		}

		crossProductGeneral(crossProductResult, vectorDifferenceA, vectorDifferenceB);

		for (int i = 0; i < DIMENSION; i++)
		{
			defVertices[tetraList[2 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			defVertices[tetraList[3 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			defVertices[tetraList[0 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[2 * numTetra + currentTetrad]] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[3 * numTetra + currentTetrad]] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[0 * numTetra + currentTetrad]] += crossProductResult[i];
		}
		tetraCounts[tetraList[2 * numTetra + currentTetrad]]++;
		tetraCounts[tetraList[3 * numTetra + currentTetrad]]++;
		tetraCounts[tetraList[0 * numTetra + currentTetrad]]++;

		//glVertex3f(defVertices[tetraList[2 * numTetra + i]].position[0], defVertices[tetraList[2 * numTetra + i]].position[1], defVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2
		//glVertex3f(defVertices[tetraList[3 * numTetra + i]].position[0], defVertices[tetraList[3 * numTetra + i]].position[1], defVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3
		//glVertex3f(defVertices[tetraList[0 * numTetra + i]].position[0], defVertices[tetraList[0 * numTetra + i]].position[1], defVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0

		for (int i = 0; i < DIMENSION; i++)
		{
			vectorDifferenceA[i] =  defVertices[tetraList[0 * numTetra + currentTetrad]].position[i] - defVertices[tetraList[1 * numTetra + currentTetrad]].position[i];
			vectorDifferenceB[i] =  defVertices[tetraList[1 * numTetra + currentTetrad]].position[i] - defVertices[tetraList[2 * numTetra + currentTetrad]].position[i];
		}

		crossProductGeneral(crossProductResult, vectorDifferenceA, vectorDifferenceB);

		for (int i = 0; i < DIMENSION; i++)
		{
			defVertices[tetraList[0 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			defVertices[tetraList[1 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			defVertices[tetraList[2 * numTetra + currentTetrad]].vertexNormal[i] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[0 * numTetra + currentTetrad]] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[1 * numTetra + currentTetrad]] += crossProductResult[i];
			//vertexNormals[i * numVertices + tetraList[2 * numTetra + currentTetrad]] += crossProductResult[i];
		}
		tetraCounts[tetraList[0 * numTetra + currentTetrad]]++;
		tetraCounts[tetraList[1 * numTetra + currentTetrad]]++;
		tetraCounts[tetraList[2 * numTetra + currentTetrad]]++;

		//glVertex3f(defVertices[tetraList[0 * numTetra + i]].position[0], defVertices[tetraList[0 * numTetra + i]].position[1], defVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0
		//glVertex3f(defVertices[tetraList[1 * numTetra + i]].position[0], defVertices[tetraList[1 * numTetra + i]].position[1], defVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
		//glVertex3f(defVertices[tetraList[2 * numTetra + i]].position[0], defVertices[tetraList[2 * numTetra + i]].position[1], defVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2


	}

	for (int i = 0; i < numVertices; i++)
	{
		double magnitude = 0;
		for (int j = 0; j < DIMENSION; j++)
		{
			defVertices[i].vertexNormal[j] /= tetraCounts[i];
			//vertexNormals[j * numVertices + i] /= tetraCounts[i];  //Find average direction of all normals for adjacent triangles
			magnitude += defVertices[i].vertexNormal[j] * defVertices[i].vertexNormal[j];
			//magnitude += vertexNormals[j * numVertices + i] * vertexNormals[j * numVertices + i];
		}

		magnitude = sqrt(magnitude);

		for (int j = 0; j < DIMENSION; j++)
		{
			defVertices[i].vertexNormal[j] /= magnitude;
			//vertexNormals[j * numVertices + i] /= magnitude;  //Normalize the normal to be unit length
		}
		
	}


//
//	int edgeCounter = 0;	//Represents current edge number
//	int i = 0, j = 0;		//For loop indexes
//	for (i = 0; i < rows - 1; i++) //row
//	{
//		for (j = 0; j < cols - 1; j++) //column
//		{
//			Particle * topLeft = &particles[edgeList[edgeCounter].start];			//because this was the start of the first connection
//			Particle * bottomLeft = &particles[edgeList[edgeCounter].end];			//because the first connection connected to the bottom
//			Particle * bottomRight = &particles[edgeList[edgeCounter + 1].end];		//because the second connection connected to the bottom right
//			Particle * topRight = &particles[edgeList[edgeCounter + 2].end];		//because the third connection connected to the right
//	
//			edgeCounter += 4;
//
//			//Corresponds to this triangle:
//			//glVertex3f(bottomLeft->position[0], bottomLeft->position[1], bottomLeft->position[2]); 
//			//glVertex3f(topLeft->position[0], topLeft->position[1], topLeft->position[2]);
//			//glVertex3f(topRight->position[0], topRight->position[1], topRight->position[2]);
//			for (int i = 0; i < DIMENSION; i++)
//			{
//				vectorDifferenceA[i] = bottomLeft -> position[i] - topLeft -> position[i];
//				vectorDifferenceB[i] = topRight -> position[i] - topLeft -> position[i];
//			}
//
//			double x = -(vectorDifferenceA[1] * vectorDifferenceB[2] - vectorDifferenceA[2] * vectorDifferenceB[1]);
//			double y = -(vectorDifferenceA[2] * vectorDifferenceB[0] - vectorDifferenceA[0] * vectorDifferenceB[2]);
//			double z = -(vectorDifferenceA[0] * vectorDifferenceB[1] - vectorDifferenceA[1] * vectorDifferenceB[0]);
//
//			bottomLeft -> normal[0] += x;
//			bottomLeft -> normal[1] += y;
//			bottomLeft -> normal[2] += z;
//			bottomLeft -> triangleCount++;
//			topLeft -> normal[0] += x;
//			topLeft -> normal[1] += y;
//			topLeft -> normal[2] += z;
//			topLeft -> triangleCount++;
//			topRight -> normal[0] += x;
//			topRight -> normal[1] += y;
//			topRight -> normal[2] += z;
//			topRight -> triangleCount++;
//					
//			//Corresponds to this triangle:
//			//glVertex3f(bottomLeft->position[0], bottomLeft->position[1], bottomLeft->position[2]); 
//			//glVertex3f(bottomRight->position[0], bottomRight->position[1], bottomRight->position[2]);
//			//glVertex3f(topRight->position[0], topRight->position[1], topRight->position[2]);
//			for (int i = 0; i < DIMENSION; i++)
//			{
//				vectorDifferenceA[i] = topRight -> position[i] - bottomRight -> position[i];
//				vectorDifferenceB[i] = bottomLeft -> position[i] - bottomRight -> position[i];
//			}
//
//			x = -(vectorDifferenceA[1] * vectorDifferenceB[2] - vectorDifferenceA[2] * vectorDifferenceB[1]);
//			y = -(vectorDifferenceA[2] * vectorDifferenceB[0] - vectorDifferenceA[0] * vectorDifferenceB[2]);
//			z = -(vectorDifferenceA[0] * vectorDifferenceB[1] - vectorDifferenceA[1] * vectorDifferenceB[0]);
//
//			bottomLeft -> normal[0] += x;
//			bottomLeft -> normal[1] += y;
//			bottomLeft -> normal[2] += z;
//			bottomLeft -> triangleCount++;
//			bottomRight -> normal[0] += x;
//			bottomRight -> normal[1] += y;
//			bottomRight -> normal[2] += z;
//			bottomRight -> triangleCount++;
//			topRight -> normal[0] += x;
//			topRight -> normal[1] += y;
//			topRight -> normal[2] += z;
//			topRight -> triangleCount++;
//			
//		}
//
//		edgeCounter++;
//	}
//
//	//For each vertex, we need a normal.  We are looking at the normal of each of the adjacent triangles and calculating an average normal.
//	//We already did the summing for this average above.  Now we just need to divide by the number of triangles.
//	//Once we find an average normal, we will then normalize it (make it unit lenth) for OpenGL's sake
//	for (int i = 0; i < numParticles; i++)
//	{
//		double magnitude = 0;
//		for (int j = 0; j < DIMENSION; j++)
//		{
//			particles[i].normal[j] /= particles[i].triangleCount;  //Find average direction of all normals for adjacent triangles
//			magnitude += particles[i].normal[j] * particles[i].normal[j];
//		}
//
//		magnitude = sqrt(magnitude);
//
//		for (int j = 0; j < DIMENSION; j++)
//		{
//			particles[i].normal[j] /= magnitude;  //Normalize the normal to be unit length
//		}
//		
//	}
//
//	#ifdef DEBUGGING
//	if (logger -> isLogging)
//	{
//		logger -> printNormals(particles,numParticles, logger -> MEDIUM);
//	}
//	#endif
}

//Render the particles
void ParticleSystem::doRender(double videoWriteDeltaT)
{
	
	//glTranslatef(-halfWidth, 0.0f, halfHeight);
	//#ifdef DEBUGGING
	//if (logger -> isLogging)
	//{
	//	cout << "Entering do Render method" << endl;
	//	logger ->printVelocitiesAndPositions(defVertices,numVertices,"Deformed Positions", "Deformed Velocities", logger ->FULL);
	///#endif

	for (int i = 0; i < numTetra; i++)
	{
		float myTetraColor(0.0);
		if (i == 0)
		{
			myTetraColor = 0.0f; // Green for tetra 1
		}
		else
		{
			//myTetraColor = 1.0f;  //Will end up with Green + Red for tetra 2
			myTetraColor = 0.0f;
		}

		if (false)
		{
			glColor3f(1.0f, myTetraColor, 0.0f);
			glBegin(GL_TRIANGLES);
			//Counter clockwise winding
			//Triangle 1
			glVertex3f(orgVertices[tetraList[3 * numTetra + i]].position[0], orgVertices[tetraList[3 * numTetra + i]].position[1], orgVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3
			glVertex3f(orgVertices[tetraList[1 * numTetra + i]].position[0], orgVertices[tetraList[1 * numTetra + i]].position[1], orgVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
			glVertex3f(orgVertices[tetraList[0 * numTetra + i]].position[0], orgVertices[tetraList[0 * numTetra + i]].position[1], orgVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0
		
			glColor3f(0.8f, myTetraColor, 0.0f);
			//Triangle 2
			glVertex3f(orgVertices[tetraList[2 * numTetra + i]].position[0], orgVertices[tetraList[2 * numTetra + i]].position[1], orgVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2
			glVertex3f(orgVertices[tetraList[1 * numTetra + i]].position[0], orgVertices[tetraList[1 * numTetra + i]].position[1], orgVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
			glVertex3f(orgVertices[tetraList[3 * numTetra + i]].position[0], orgVertices[tetraList[3 * numTetra + i]].position[1], orgVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3

			glColor3f(0.6f, myTetraColor, 0.0f);
			//Triangle 3
			glVertex3f(orgVertices[tetraList[2 * numTetra + i]].position[0], orgVertices[tetraList[2 * numTetra + i]].position[1], orgVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2
			glVertex3f(orgVertices[tetraList[3 * numTetra + i]].position[0], orgVertices[tetraList[3 * numTetra + i]].position[1], orgVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3
			glVertex3f(orgVertices[tetraList[0 * numTetra + i]].position[0], orgVertices[tetraList[0 * numTetra + i]].position[1], orgVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0

			glColor3f(0.4f, myTetraColor, 0.0f);
			//Triangle 4
			glVertex3f(orgVertices[tetraList[0 * numTetra + i]].position[0], orgVertices[tetraList[0 * numTetra + i]].position[1], orgVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0
			glVertex3f(orgVertices[tetraList[1 * numTetra + i]].position[0], orgVertices[tetraList[1 * numTetra + i]].position[1], orgVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
			glVertex3f(orgVertices[tetraList[2 * numTetra + i]].position[0], orgVertices[tetraList[2 * numTetra + i]].position[1], orgVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2

			glEnd();
		}

		glColor3f(myTetraColor, 1.0, 0.0f);
		glBegin(GL_TRIANGLES);
		
		
		int normalVertex = 0;
		//Counter clockwise winding
		//Triangle 1

		normalVertex = tetraList[3 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[3 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[3 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[3 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[3 * numTetra + i]].position[0], defVertices[tetraList[3 * numTetra + i]].position[1], defVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3
		normalVertex = tetraList[1 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[1 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[1 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[1 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[1 * numTetra + i]].position[0], defVertices[tetraList[1 * numTetra + i]].position[1], defVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
		normalVertex = tetraList[0 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[0 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[0 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[0 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[0 * numTetra + i]].position[0], defVertices[tetraList[0 * numTetra + i]].position[1], defVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0
		

		glColor3f(myTetraColor, 0.8f, 0.0f);
		normalVertex = tetraList[2 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[2 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[2 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[2 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[2 * numTetra + i]].position[0], defVertices[tetraList[2 * numTetra + i]].position[1], defVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2
		normalVertex = tetraList[1 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[1 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[1 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[1 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[1 * numTetra + i]].position[0], defVertices[tetraList[1 * numTetra + i]].position[1], defVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
		normalVertex = tetraList[3 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[3 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[3 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[3 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[3 * numTetra + i]].position[0], defVertices[tetraList[3 * numTetra + i]].position[1], defVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3


		glColor3f(myTetraColor, 0.6f, 0.0f);
		normalVertex = tetraList[2 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[2 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[2 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[2 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[2 * numTetra + i]].position[0], defVertices[tetraList[2 * numTetra + i]].position[1], defVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2
		normalVertex = tetraList[3 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[3 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[3 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[3 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[3 * numTetra + i]].position[0], defVertices[tetraList[3 * numTetra + i]].position[1], defVertices[tetraList[3 * numTetra + i]].position[2]);  //vertex 3
		normalVertex = tetraList[0 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[0 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[0 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[0 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[0 * numTetra + i]].position[0], defVertices[tetraList[0 * numTetra + i]].position[1], defVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0


		glColor3f(myTetraColor, 0.4f, 0.0f);
		normalVertex = tetraList[0 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[0 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[0 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[0 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[0 * numTetra + i]].position[0], defVertices[tetraList[0 * numTetra + i]].position[1], defVertices[tetraList[0 * numTetra + i]].position[2]);  //vertex 0
		normalVertex = tetraList[1 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[1 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[1 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[1 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[1 * numTetra + i]].position[0], defVertices[tetraList[1 * numTetra + i]].position[1], defVertices[tetraList[1 * numTetra + i]].position[2]);  //vertex 1
		normalVertex = tetraList[2 * numTetra + i];
		glNormal3f(defVertices[normalVertex].vertexNormal[0], defVertices[normalVertex].vertexNormal[1], defVertices[normalVertex].vertexNormal[2]);
		//glNormal3f(vertexNormals[numVertices * 0 + tetraList[2 * numTetra + i]], vertexNormals[numVertices * 1 + tetraList[2 * numTetra + i]], vertexNormals[numVertices * 2 + tetraList[2 * numTetra + i]]);
		glVertex3f(defVertices[tetraList[2 * numTetra + i]].position[0], defVertices[tetraList[2 * numTetra + i]].position[1], defVertices[tetraList[2 * numTetra + i]].position[2]);  //vertex 2

		glEnd();
		
	}

	//Draw Normals
	if (showNormals)
	{
		int verticesForNormals [3][4] = {{0, 1, 0 , 0}, {1, 2, 1, 2}, {2, 3, 3, 3}};
		double midPoint[3];
		double midPlusNormal[3];
		for (int i = 0; i < numTetra; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				midPoint[0] = (defVertices[tetraList[verticesForNormals[0][j] * numTetra+ i]].position[0] + defVertices[tetraList[verticesForNormals[1][j] * numTetra+ i]].position[0] + defVertices[tetraList[verticesForNormals[2][j] * numTetra+ i]].position[0]) / 3.0;
				midPoint[1] = (defVertices[tetraList[verticesForNormals[0][j] * numTetra+ i]].position[1] + defVertices[tetraList[verticesForNormals[1][j] * numTetra+ i]].position[1] + defVertices[tetraList[verticesForNormals[2][j] * numTetra+ i]].position[1]) / 3.0;
				midPoint[2] = (defVertices[tetraList[verticesForNormals[0][j] * numTetra+ i]].position[2] + defVertices[tetraList[verticesForNormals[1][j] * numTetra+ i]].position[2] + defVertices[tetraList[verticesForNormals[2][j] * numTetra+ i]].position[2]) / 3.0;

				//Vertex vertex1 = defVertices[tetraList[verticesForNormals[0][j] * numTetra+ i]];
				//cout << "Midpoint vertex #1:" << endl;
				//for (int i = 0; i < 3; i++)
				//{
				//	cout << vertex1.position[i] << " ";
				//}
				//cout << endl;

				//Vertex vertex2 = defVertices[tetraList[verticesForNormals[1][j] * numTetra+ i]];
				//cout << "Midpoint vertex #2:" << endl;
				//for (int i = 0; i < 3; i++)
				//{
				//	cout << vertex2.position[i] << " ";
				//}
				//cout << endl;
				
				//Vertex vertex3 = defVertices[tetraList[verticesForNormals[2][j] * numTetra+ i]];
				//cout << "Midpoint vertex #3:" << endl;
				//for (int i = 0; i < 3; i++)
				//{
				//	cout << vertex3.position[i] << " ";
				//}
				//cout << endl;

				#ifdef DEBUGGING
				if (logger -> isLogging)
				{
					logger -> printVector(midPoint,3,"midPoint", logger -> FULL);
					cout << endl;
				}
				#endif


				midPlusNormal[0] = midPoint[0] + normals[0 * 4 * numTetra + i * 4 + j];
				midPlusNormal[1] = midPoint[1] + normals[1 * 4 * numTetra + i * 4 + j];
				midPlusNormal[2] = midPoint[2] + normals[2 * 4 * numTetra + i * 4 + j];

				#ifdef DEBUGGING
				if (logger -> isLogging)
				{
					logger -> printVector(midPlusNormal,3,"midPlusNormal", logger -> FULL);
					cout << endl;
				}
				#endif

				glColor3f(1.0f, 0.0f, 0.0f);
				glBegin(GL_LINES);
				glVertex3f(midPoint[0], midPoint[1], midPoint[2]);
				glVertex3f(midPlusNormal[0], midPlusNormal[1], midPlusNormal[2]);
				glEnd();
			}
		}
	}

	glBegin(GL_QUADS);
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(-10.0f, -4.0f, -10.0f);
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(10.0f, -4.0f, -10.0f);
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(10.0f, -4.0f, 10.0f);
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f(-10.0f, -4.0f, 10.0f);
	glEnd();

	//Render onscreen text with informational messages
	//Note: This now will be written to the video images
	if (showInfoText)
	{
		//Reference: http://programming-technique.blogspot.com/2012/05/font-rendering-in-glut-using-bitmap.html
		glDisable(GL_LIGHTING);
		glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, windowWidth, 0, windowHeight);
		
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glRasterPos2i(20, 20);
		for (int i = 0; i < strlen(text); i++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, text[i]);
		}
		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		glEnable(GL_LIGHTING);
	}

	if (renderToImage)
	{

		if (timeSinceVideoWrite >= videoWriteDeltaT)
		{
			timeSinceVideoWrite = 0.0;
			//Read the rendered image into a buffer
			glFlush();
			glReadPixels(0, 0, windowWidth, windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, screenBuffer);

			//Vertically flip the image in memory - it flips two rows from the screen in memory at a time
			for (int i = 0; i < windowHeight / 2; i++)
			{
				if (i != windowHeight - i)
				{
					memcpy(screenRowTemp,(uint8_t *)&screenBuffer[i * windowWidth * 4], windowWidth * 4 * sizeof(uint8_t));
					memcpy((uint8_t *)&screenBuffer[i * windowWidth * 4], (uint8_t *)&screenBuffer[(windowHeight - i) * windowWidth * 4], windowWidth * 4 * sizeof(uint8_t));
					memcpy((uint8_t *)&screenBuffer[(windowHeight - i) * windowWidth * 4], screenRowTemp, windowWidth * 4 * sizeof(uint8_t));
				}
			}
		
			//Write a new targa file, with a new number
			sprintf(imageFileName, "images\\ImplicitMethods%d.tga", frameNumber);
			tga_result result = tga_write_rgb(imageFileName, screenBuffer, windowWidth, windowHeight, 32);
			#ifdef DEBUGGING
			if (logger -> isLogging && logger -> loggingLevel >= logger -> LIGHT)
			{
				if (result != TGA_NOERR)
				{
					cout << tga_error(result) << " at frame number " << frameNumber << endl;
				}
			}
			#endif

			frameNumber++;
		}

		
	}
		
}

//This method stores the current eye position (from the camera) for this mesh
//parameter eyePos - the eye position to store
void ParticleSystem::setEyePos(glm::vec3 & eyePos)
{
	this -> eyePos[0] = eyePos[0];
	this -> eyePos[1] = eyePos[1];
	this -> eyePos[2] = eyePos[2];
}

void ParticleSystem::initVBOs()
{
	glGenBuffers(1, vboHandle);
	glGenBuffers(1, indexVboHandle);

	glGenBuffers(1, floorVboHandle);
	glGenBuffers(1, floorIndexVboHandle);
}

void ParticleSystem::sendVBOs()
{
	//Tetrahedral Mesh
	glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * numVertices, defVertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	//glColor3f(myTetraColor, 1.0, 0.0f);
	indices.clear();
	for (int i = 0; i < numTetra; i++)
	{
		//Counter clockwise winding
		//Triangle 1
		indices.push_back(tetraList[3 * numTetra + i]);
		indices.push_back(tetraList[1 * numTetra + i]);
		indices.push_back(tetraList[0 * numTetra + i]);
		
		indices.push_back(tetraList[2 * numTetra + i]);
		indices.push_back(tetraList[1 * numTetra + i]);
		indices.push_back(tetraList[3 * numTetra + i]);

		indices.push_back(tetraList[2 * numTetra + i]);
		indices.push_back(tetraList[3 * numTetra + i]);
		indices.push_back(tetraList[0 * numTetra + i]);

		indices.push_back(tetraList[0 * numTetra + i]);
		indices.push_back(tetraList[1 * numTetra + i]);
		indices.push_back(tetraList[2 * numTetra + i]);

	}

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVboHandle[0]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * indices.size(), &indices[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	//Floor
	floorVertices.clear();
	//Vertex vertex0 = {-halfWidth, -halfHeight, halfDepth, 1, normal0[0], normal0[1], normal0[2], 0, color1[0], color1[1], color1[2], 1.0f}; //Front lower left
	Vertex vertex0 = Vertex(-10.0f, -4.0f, -10.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0);
	Vertex vertex1 = Vertex(10.0f, -4.0f, -10.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f,  0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0);
	Vertex vertex2 = Vertex(10.0f, -4.0f, 10.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0);
	Vertex vertex3 = Vertex(-10.0f, -4.0f, 10.0f, 1.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0);

	floorVertices.push_back(vertex0);
	floorVertices.push_back(vertex1);
	floorVertices.push_back(vertex2);
	floorVertices.push_back(vertex3);

	//glBegin(GL_QUADS);
		//glNormal3f(0.0f, 1.0f, 0.0f);
		//glVertex3f(-10.0f, -4.0f, -10.0f);
		//glNormal3f(0.0f, 1.0f, 0.0f);
		//glVertex3f(10.0f, -4.0f, -10.0f);
		//glNormal3f(0.0f, 1.0f, 0.0f);
		//glVertex3f(10.0f, -4.0f, 10.0f);
		//glNormal3f(0.0f, 1.0f, 0.0f);
		//glVertex3f(-10.0f, -4.0f, 10.0f);
	//glEnd();

	//Data structure:
	//float position [DIMENSION + 1];	//Position vector in 3d space
	//float vertexNormal [DIMENSION + 1];		//Normal vector in 3d space (this one is a VERTEX AVERAGE used in lighting)
	//float color[4];
	//float velocity [DIMENSION + 1];	//Velocity vector in 3d space	
	//int triangleCount;				//Number of triangles in which the vertex for this particle is contained (helps calculate normals)
	
	floorIndices.clear();
	//Two triangles for a quad
	floorIndices.push_back(0);
	floorIndices.push_back(1);
	floorIndices.push_back(2);
	floorIndices.push_back(2);
	floorIndices.push_back(3);
	floorIndices.push_back(0);

	glBindBuffer(GL_ARRAY_BUFFER, floorVboHandle[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex) * floorVertices.size(), &floorVertices[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, floorIndexVboHandle[0]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * floorIndices.size(), &floorIndices[0], GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);


}

void ParticleSystem::doRender(double videoWriteDeltaT, glm::mat4 & projMatrix, glm::mat4 & modelViewMatrix)
{
	sendVBOs();

	glm::mat4 totalMatrix = projMatrix * modelViewMatrix;
	glm::mat4 normalMatrix = glm::inverse(modelViewMatrix);
	normalMatrix = glm::transpose(normalMatrix);

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

	glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(Vertex),(char*) NULL+0); 
	glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(Vertex),(char*) NULL+16); 
    glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(Vertex),(char*) NULL+32); 

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


	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////
	{

	glm::mat4 totalMatrix = projMatrix * modelViewMatrix;
	glm::mat4 normalMatrix = glm::inverse(modelViewMatrix);
	normalMatrix = glm::transpose(normalMatrix);

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

	glBindBuffer(GL_ARRAY_BUFFER, floorVboHandle[0]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, floorIndexVboHandle[0]);

	glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(Vertex),(char*) NULL+0); 
	glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(Vertex),(char*) NULL+16); 
    glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(Vertex),(char*) NULL+32); 

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
	glDrawElements(GL_TRIANGLES, floorIndices.size(), GL_UNSIGNED_INT, (char *) NULL + 0);

	glUseProgram(0);
	

	}

	

	
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////User Interface methods to alter attributes of the particle system///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Method to increase the earth gravity used in the particle system
//Parameter amount - the amount to increase.  Use a negative value to decrease
void ParticleSystem::increaseEarthGravity(double amount)
{
	earthGravityValue += amount;
	if (true)
	{
		cout << "Earth Gravity value is now " << earthGravityValue << endl;
	}

	sprintf(text, "Earth Gravity: %f", earthGravityValue);
}

//Method to increase the straight rest length
//Parameter amount - the amount by which to increase it (a negative value decreases)
void ParticleSystem::increaseStraightRestLength(double amount)
{
	//restLengthValues[0] += amount;
	//restLengthValues[1] += amount;
	//if (restLengthValues[0] < epsilon || restLengthValues[1] < epsilon)	//Rest length of 0 results in bad calculation of the b vector.  Negative value doesn't make sense.
 //   {
	//	restLengthValues[0] -= amount;
	//	restLengthValues[1] -= amount;
	//	return;
 //   }

	//restLengthValues[2] = sqrt(restLengthValues[0] * restLengthValues[0] + restLengthValues[1] * restLengthValues[1]);  //Use pythagorean theorem to calculate diagonal restLength
	//

	////Half of the width and half of the height have now changed and must be calculated again
	//sheetWidth += amount * (cols - 1);
	//sheetHeight += amount * (rows - 1);
	//halfWidth = sheetWidth / 2;
	//halfHeight = sheetHeight / 2;

	////Recompute all constrained particles in the grid.
	////Non constrained particles will be allowed to "fall in line" due to their spring forces
	//		
	//int currentParticle = 0;
	//for (int i = 0; i < rows; i++)
	//{
	//	for (int j = 0; j < cols; j++)
	//	{
	//		currentParticle = i * cols + j;
	//		for (int k = 0; k < numConstraints; k++)
	//		{
	//			if (currentParticle == constraintParticles[k])
	//			{
	//				particles[currentParticle].position[0] = 0 + j * restLengthValues[0];
	//				particles[currentParticle].position[2] = 0 + -i * restLengthValues[1];
	//				particles[currentParticle].position[1] = 0;
	//				//particles[i * cols + j].position.y = 0 + -i * padding;		//Useful for a vertical sheet
	//				//particles[i * cols + j].position.z = 0;						//Useful for a vertical sheet
	//				for (int m = 0; m < DIMENSION; m++)
	//				{
	//					particles[currentParticle].velocity[m] = 0;
	//				}
	//				break;
	//			}
	//		}

	//	}
	//}

	//if (true)
	//{
	//	cout << "Straight Rest Length is now  " << restLengthValues[0] << endl;
	//	cout << "Diagonal Rest Length is now " << restLengthValues[1] << endl;
	//}

	//sprintf(text, "Straight Rest Length: %f", restLengthValues[0]);
}

//Method to toggle whether or not informational messages are shown on screen
void ParticleSystem::toggleInfoText()
{
	showInfoText = !showInfoText;
}

//Method to set the window width and height
//Parameters width and height - the window dimensions
void ParticleSystem::setWindowDimensions(int width, int height)
{
	windowWidth = width;
	windowHeight = height;
}

//Method to toggle whether or not animation occurs
//Particles will not have their positions or velocities updated if animation is turned off
void ParticleSystem::toggleAnimation()
{
	isAnimating = !isAnimating;
}

//Method to toggle between rendering surfaces and wireframe mode
void ParticleSystem::toggleRenderMode()
{
	renderMode = !renderMode;
}

//Method to toggle between rendering to a series of numbered images and not rendering to them
//It ALWAYS renders to the screen regardless
void ParticleSystem::toggleImageRendering()
{
	renderToImage = !renderToImage;
	if (renderToImage)
	{
		//sprintf(text, "Image rendering on");
	}
	else
	{
		//sprintf(text, "Image rendering off");
	}
}

void ParticleSystem::doTransform()
{

	//return;
	
	double identity[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
	double scaleFactor = 1;
	double scaleTransform[3][3] = {{scaleFactor, 0, 0}, {0, scaleFactor, 0}, {0, 0, scaleFactor}};
	double angle = -90;
	//double angle = 0;
	double rotateTransform[3][3] = {{1, 0, 0}, {0, cos (angle * 3.14 / 180), -sin (angle * 3.14 / 180)}, {0, sin (angle * 3.14 / 180), cos (angle * 3.14 / 180)}};
	double totalTransform[3][3];
	
	multiplyMatricesATimesB(totalTransform,rotateTransform,3,3,scaleTransform,3,3);  //couldn't use because it requires two subscripts
	
	//multiplyMatricesATimesB(totalTransform,scaleTransform,3,3,rotateTransform,3,3);  //couldn't use because it requires two subscripts

	////////////////////////////////
	//do transform

	
	double * temp = new double [DIMENSION * numVertices];
	
	for (int i = 0; i < DIMENSION; i++) //Current row # of calculated final matrix element						
	{																											
		for (int j = 0; j < numVertices; j++) //current col # of calculated final matrix element					
		{																										
			temp[i * numVertices + j] = 0;																				
			for (int k = 0; k < 3; k++) //current col in left matrix and row in right matrix for summing
			{																									
				temp[i * numVertices + j] += totalTransform[i][k] * defVertices[j].position[k];											
			}																									
																												
		}																										
																												
	}			
	
	for (int i = 0; i < DIMENSION; i++)
	{
		for (int j = 0; j < numVertices; j++)
		{
			defVertices[j].position[i] = temp[i * numVertices + j];
		}
	}

	delete [] temp;

}

void ParticleSystem::printStateReport()
{
	/*
	for (int i = 0; i < numVertices; i++)
	{
		cout << "Vertex " << i << " Position:" << endl;
		cout << defVertices[i].position[0] << " " << defVertices[i].position[1] << " " << defVertices[i].position[2] << endl;
	}

	for (int i = 0; i < numVertices; i++)
	{
		cout << "Vertex " << i << " Velocity :" << endl;
		cout << defVertices[i].velocity[0] << " " << defVertices[i].velocity[1] << " " << defVertices[i].velocity[2] << endl;
	}
	*/

	cout << "lambda: " << lambda << endl;
	cout << "mu: " << mu << endl;
	cout << "kd: " << kd << endl;

	/*
	for (int i = 0; i < numVertices; i++)
	{
		cout << "Force for vertex " << i << ":" << endl;
		for (int j = 0; j < DIMENSION; j++)
		{
			cout << currentForce[i * DIMENSION + j] << " ";
		}
		cout << endl;
	}*/

	for (int i = 0; i < numTetra; i++)
	{
		cout << "Tetra " << i << ":" << endl;
		for (int j = 0; j < 4; j++)
		{
			cout << "Vertex " << j << ":" << endl;
			///cout << defVertices[tetraList[
			cout << defVertices[tetraList[j * numTetra + i]].position[0] << " " << defVertices[tetraList[j * numTetra + i]].position[1] << " " << defVertices[tetraList[j * numTetra + i]].position[2] << endl;;
		}
	}


}

//This method is going to load up a special state that lets us get roll'n with some error anlaysis!
void ParticleSystem::loadSpecialState()
{
	
	//Vertex 0 Position:
	//-0.00664436 1.4789 -0.777403
	//Vertex 1 Position:
	//2.15984 -0.137176 -2.39931
	//Vertex 2 Position:
	//1.00832 -0.99262 0.434279
	//Vertex 3 Position:
	//-0.472082 -0.970312 -2.10122
	//Vertex 4 Position:
	//2.79938 -0.999992 -1.83487
	//Vertex 0 Velocity :
	//0.476164 0.18778 -0.121715
	//Vertex 1 Velocity :
	//-0.0326747 -1.39825 -0.148798
	//Vertex 2 Velocity :
	//-0.052081 0.126357 0.102539
	//Vertex 3 Velocity :
	//-0.132509 0.305334 0.088056
	//Vertex 4 Velocity :
	//0.00300343 0.015784 -0.00969502
	//lambda: 116.667
	//mu: 350
	//kd: 10
	//Force for vertex 0:
	//0.737907 -4.87486 1.83629
	//Force for vertex 1:
	//1.49416 0.0737536 16.8846
	//Force for vertex 2:
	//6.60933 9.91677 9.63472
	//Force for vertex 3:
	//-21.9418 1.47707 2.00352
	//Force for vertex 4:
	//-1.07333 -1.2923 -1.33103
	
	

	
	  

	defVertices[0].position[0] = -0.00664436;
	defVertices[0].position[1] = 1.4789;
	defVertices[0].position[2] = -0.777403;
	defVertices[1].position[0] = 2.15984;
	defVertices[1].position[1] = -0.137176;
	defVertices[1].position[2] = -2.39931;
	defVertices[2].position[0] = 1.00832;
	defVertices[2].position[1] = -0.99262;
	defVertices[2].position[2] = 0.434279;
    defVertices[3].position[0] = -0.472082;
	defVertices[3].position[1] = -0.970312;
	defVertices[3].position[2] = -2.10122;
	defVertices[4].position[0] = 2.79938;
	defVertices[4].position[1] = -0.999992;
	defVertices[4].position[2] = -1.83487;

	

	defVertices[0].velocity[0] = 0.476164;
	defVertices[0].velocity[1] = 0.18778;
	defVertices[0].velocity[2] = -0.121715;
	defVertices[1].velocity[0] = -0.0326747;
	defVertices[1].velocity[1] = -1.39825;
	defVertices[1].velocity[2] = -0.148798;
	defVertices[2].velocity[0] = -0.052081;
	defVertices[2].velocity[1] = 0.126357;
	defVertices[2].velocity[2] = 0.102539;
    defVertices[3].velocity[0] = -0.132509;
	defVertices[3].velocity[1] = 0.305334;
	defVertices[3].velocity[2] = 0.088056;
	defVertices[4].velocity[0] = 0.00300343;
	defVertices[4].velocity[1] = 0.015784;
	defVertices[4].velocity[2] = -0.00969502;

	lambda = 116.667;
	mu = 350;
	kd = 10;

	


}