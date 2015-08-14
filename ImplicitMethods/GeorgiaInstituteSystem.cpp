#include <iostream>
#include <iomanip>
#include <assert.h>
#include "Logger.h"
#include "GeorgiaInstituteSystem.h"

using namespace std;

#include "Macros.h"

const double epsilon = 1e-12;	//Used to check approximate equality to 0

//Based on the paper at http://graphics.berkeley.edu/papers/Obrien-GMA-1999-08/Obrien-GMA-1999-08.pdf – Graphical Modeling and Animation of Brittle Fracture
//By James O'Brien and Jessica Hodgkins

//Constructor
GeorgiaInstituteSystem::GeorgiaInstituteSystem(Vertex * vertexList, int vertexCount, int * tetraList, int tetraCount, Logger * logger) : ParticleSystem(vertexList, vertexCount, tetraList, tetraCount, logger)
{
	strcpy(text, "Method 2");
	////double K = 100;					//Bulk Modulus
	////mu = 100;						//Shear modulus (Lame's second parameter)
	double K = 700;
	mu = 700;
	//double K = 20;								//Bulk Modulus
	//mu = 20;	
	lambda = K - (2.0/3) * mu;		//Lame's first parameter
	//kd = 1.93;
	kd = 0.2;

	//doTransform();

	cout << fixed << setprecision(5);
	logger -> initFile("georgiaDeformation.log");
	logger ->printText(text);

	//m = zeros(3,size(triangles,2)*4);
	m = new double[4 * 4 * numTetra];
	//beta = zeros(4,size(triangles,2)*4);
	beta = new double [4 * 4 * numTetra];

	for (int currentTetrad = 0; currentTetrad < numTetra; currentTetrad++)
	{
		for(int j = 0; j < DIMENSION; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				m[4 * numTetra * j + currentTetrad * 4 + k] = orgVertices[tetraList[k * numTetra + currentTetrad]].position[j];
				
			}
			
		}

		for (int j = 0; j < 4; j++)
		{
			m[12 * numTetra + currentTetrad * 4 + j] = 1;
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print4By4MatrixSingleIndex(m, numTetra, currentTetrad, "m", logger -> MEDIUM);
		}
		#endif

		//Formula from http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
		double invDeterminant =
		1.0 /
		(
		  m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+1] + m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+2]
		+ m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+2] + m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+0]
		+ m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+0] + m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+1]
		+ m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+1] + m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+2] + m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+0]
		- m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+1]
		- m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+2]
		- m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+0]
		- m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+1]
		);

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger ->loggingLevel >= logger->MEDIUM)
			{
				logger->printIteration("Inverse of determinant is: ",invDeterminant);
			}
		}
		#endif

		beta[0*numTetra*4+currentTetrad*4+0] = invDeterminant * (m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+3] + m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+1] + m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+2] - m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+2] - m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+3] - m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+1]);
		beta[0*numTetra*4+currentTetrad*4+1] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+2] + m[0*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+2]);
		beta[0*numTetra*4+currentTetrad*4+2] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+1] + m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+1]);
		beta[0*numTetra*4+currentTetrad*4+3] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+2] + m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+2]);
		beta[1*numTetra*4+currentTetrad*4+0] = invDeterminant * (m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+2] + m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+3] + m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+0] - m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+3] - m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+0] - m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+2]);
		beta[1*numTetra*4+currentTetrad*4+1] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+0] + m[0*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+0]);
		beta[1*numTetra*4+currentTetrad*4+2] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+2] + m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+2]);
		beta[1*numTetra*4+currentTetrad*4+3] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+0] + m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+0]);
		beta[2*numTetra*4+currentTetrad*4+0] = invDeterminant * (m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+3] + m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+0] + m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+1] - m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+1] - m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+3] - m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+0]);
		beta[2*numTetra*4+currentTetrad*4+1] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+1] + m[0*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+1]);
		beta[2*numTetra*4+currentTetrad*4+2] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+0] + m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+3]*m[3*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+0]);
		beta[2*numTetra*4+currentTetrad*4+3] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+1] + m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+3] + m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+3] - m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+3]*m[2*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+3]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+1]);
		beta[3*numTetra*4+currentTetrad*4+0] = invDeterminant * (m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+1] + m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+2] + m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+0] - m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+2] - m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+0] - m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+1]);
		beta[3*numTetra*4+currentTetrad*4+1] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+2] + m[0*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+0] + m[0*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+0]);
		beta[3*numTetra*4+currentTetrad*4+2] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+1] + m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+2] + m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+1]*m[3*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+2]*m[3*numTetra*4+currentTetrad*4+0] - m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+0]*m[3*numTetra*4+currentTetrad*4+1]);
		beta[3*numTetra*4+currentTetrad*4+3] = invDeterminant * (m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+2] + m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+0] + m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+0]*m[1*numTetra*4+currentTetrad*4+2]*m[2*numTetra*4+currentTetrad*4+1] - m[0*numTetra*4+currentTetrad*4+1]*m[1*numTetra*4+currentTetrad*4+0]*m[2*numTetra*4+currentTetrad*4+2] - m[0*numTetra*4+currentTetrad*4+2]*m[1*numTetra*4+currentTetrad*4+1]*m[2*numTetra*4+currentTetrad*4+0]);

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print4By4MatrixSingleIndex(beta, numTetra, currentTetrad, "beta", logger -> MEDIUM);
		}
		#endif



	}
	
}

GeorgiaInstituteSystem::~GeorgiaInstituteSystem()
{
	delete [] m;
	delete [] beta;
	

}

int iter = 1;

//Overriden doUpdate method
void GeorgiaInstituteSystem::doUpdate(double deltaT)
{
	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printText("=======================================================================================================");
		logger -> printIteration("Iteration #", iter);
		logger ->printVelocitiesAndPositions(defVertices, numVertices, "Positions - start of method", "Velocities - start of method", logger ->MEDIUM);
	}
	#endif
	

	//Start with 0 force each iteration
	//currentForce = zeros(3, size(inVelocities,2));
	for (int i = 0; i < DIMENSION * numVertices; i++)
	{
		currentForce[i] = 0;
	}
	
	//STart
	//for i = 1:size(triangles,2)
	for (int currentTetrad = 0; currentTetrad < numTetra; currentTetrad++)
	{
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> printText("**************");
			logger -> printIteration("Tetrahedron #", currentTetrad);
		}
		#endif
		
		//3 X 4 matrix followed by an optional row of all one's (4 X 4 total)
		//double m [16];
		//3 X 4 matrices
		double p [12];
		double v [12];
		
		for(int j = 0; j < DIMENSION; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				p[j * 4 + k] = defVertices[tetraList[k * numTetra + currentTetrad]].position[j];
				v[j * 4 + k] = defVertices[tetraList[k * numTetra + currentTetrad]].velocity[j];
			}
		}

		//m = [orgVertices(:, triangles(1, i)) orgVertices(:, triangles(2, i)) orgVertices(:, triangles(3, i)) orgVertices(:, triangles(4, i))];
		//p = [defVertices(:, triangles(1, i)) defVertices(:, triangles(2, i)) defVertices(:, triangles(3, i)) defVertices(:, triangles(4, i))];
		//v = [inVelocities(:, triangles(1, i)) inVelocities(:, triangles(2, i)) inVelocities(:, triangles(3, i)) inVelocities(:, triangles(4, i))];
    
		//if isDebugging == 1
		//	display('m:');
		//	display(m);
		//	display('p:');
		//	display(p);
		//	display('v:');
		//	display(v);
		//end

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print3By4MatrixSingleIndex(p,"p", logger ->MEDIUM);
			logger -> print3By4MatrixSingleIndex(v,"v", logger ->MEDIUM);
		}
		#endif


		//beta = inv([m; 1 1 1 1]);
		
		//double beta[16];

		
		
		//if isDebugging == 1
		//	display('beta:');
		//	display(beta);
		//end
    

		//////////////////////////////
		//////////////////////////////
		//p * beta
		//3 X 4(p) * 4 X 4(beta) = 3 X 4(result) size matrix

		double resultMatrix [12];
																																									
		//for (int i = 0; i < 3; i++) //current row # of calculated final matrix element						
		//{																											
		//	for (int j = 0; j < 4; j++) //current col # of calculated final matrix element					
		//	{																										
		//		resultMatrix[i * 4 + j] = 0;																				
		//		for (int k = 0; k < 4; k++) //current col in left matrix and row in right matrix for summing
		//		{																									
		//			resultMatrix[i * 4 + j] += p[i * 4 + k] * beta[k * 4 + j];											
		//		}																									
		//																										
		//	}																										
		//																										
		//}

		//multiplyMatricesATimesBSingleIndexes(resultMatrix, p, 3, 4, beta, 4, 4);
		
		//TO DO (TODO): use a macro here that works with an offset
		int aRows = 3, aCols = 4, bCols = 4;
		
		/////////////
		for (int i = 0; i < aRows; i++) /*current row # of calculated final matrix element*/						
		{																											
			for (int j = 0; j < bCols; j++) /*current col # of calculated final matrix element*/					
			{																										
				resultMatrix[i * bCols + j] = 0;																				
				for (int k = 0; k < aCols; k++) /*current col in left matrix and row in right matrix for summing*/	
				{																									
					resultMatrix[i * bCols + j] += p[i * aCols + k] * beta[k * 4 * numTetra + currentTetrad * 4 + j];							
				}																									
																												
			}																										
																												
		}		
		/////////////

		double partialXWrtUi [3];
		double partialXWrtUj [3];
		double fullPartialXWrtU[3][3];

		

		//I = eye(3);
		//e = zeros(3,3);
		double e[9]; //3 X 3 matrix
		for (int ii = 0; ii < 3; ii++)
		{
		//for ii = 1:3
			//partialXWrtUi = p * beta * [I(ii,:) 0]';
			for (int i = 0; i < 3; i++)
			{
				partialXWrtUi[i] = 0;
				for (int j = 0; j < 4; j++)
				{
					partialXWrtUi[i] += resultMatrix[i * 4 + j] * (ii == j);
				}
				fullPartialXWrtU[ii][i] = partialXWrtUi[i]; //is it really horizontal?  Or is it vertical?
			}

			
			#ifdef DEBUGGING
			if (logger -> isLogging)
			{
				logger -> printVector(partialXWrtUi,3,"partialXWrtUi", logger ->FULL);
			}
			#endif


			//for jj = 1:3
			for (int jj = 0; jj < 3; jj++)
			{
				//partialXWrtUj = p * beta * [I(jj,:) 0]';
				for (int i = 0; i < 3; i++)
				{
					partialXWrtUj[i] = 0;
					for (int j = 0; j < 4; j++)
					{
						partialXWrtUj[i] += resultMatrix[i * 4 + j] * (jj == j);
					}
				}

				#ifdef DEBUGGING
				if (logger -> isLogging)
				{
					logger -> printVector(partialXWrtUj,3,"partialXWrtUj", logger ->FULL);
				}
				#endif

				//e(ii,jj) = dot(partialXWrtUi, partialXWrtUj) - I(ii,jj);
	            e[ii * 3 + jj] = dot3(partialXWrtUi, partialXWrtUj) - (ii == jj);
		    //end
			}
		//end
		}

		//Inverted tetrahedron code is likely to start here...
		//double thatDeterminant = determinant3By3(fullPartialXWrtU);
		
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			//setprecision(18);
			logger -> print3By3MatrixSingleIndex(e,"e", logger ->MEDIUM);
		}
		#endif

		
		//elasticStress = lambda * trace(e) * I + 2 * mu * e;

		double elasticStress[9];

		double trace = trace3By3(e);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				elasticStress[i*3+j] = lambda*trace*(i==j) + 2 * mu * e[i*3+j];

			}
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			cout << setprecision(18);
			logger -> print3By3MatrixSingleIndex(elasticStress,"elasticStress", logger ->MEDIUM);
		}
		#endif

		//volume = (1/6) * dot(crossProduct((m(:,2) - m(:,1)),(m(:,3) - m(:,1))), (m(:,4) - m(:,1)));

		//logger -> print3By4MatrixSingleIndex(m,"m", logger ->MEDIUM);

		double temp1[3];
		double temp2[3];
		double temp3[3];
		
		for (int i = 0; i < 3; i++)
		{
			temp1[i] = m[i*numTetra*4+currentTetrad*4+1] - m[i*numTetra*4+currentTetrad*4+0];
			temp2[i] = m[i*numTetra*4+currentTetrad*4+2] - m[i*numTetra*4+currentTetrad*4+0];
			temp3[i] = m[i*numTetra*4+currentTetrad*4+3] - m[i*numTetra*4+currentTetrad*4+0];
		}

		double temp4[3];
		///crossProductGeneral(temp4,temp1,temp2);
		crossProductGeneral(temp4,temp2,temp1);
		
		//logger ->printVector(temp1, 3, "temp1", logger ->LIGHT);
		//logger ->printVector(temp2, 3, "temp2", logger ->LIGHT);
		//logger ->printVector(temp3, 3, "Ttemp3", logger ->LIGHT);
		//logger ->printVector(temp4, 3, "That cross product", logger ->LIGHT);

		double volume = (1.0/6) * dot3(temp4,temp3);

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger ->loggingLevel >= logger->MEDIUM)
			{
				logger ->printIteration("Volume is:", volume);
			}
		}
		#endif
		    
		//forces = zeros(3,4);
		//for ii = 1:4
		//	%forces(:,ii) = zeros(3,1);
		//	for j = 1:4
		//		sum = 0;
		//		for k = 1:3
		//			for l = 1:3
		//				sum = sum + beta(j,l) * beta(ii,k) * elasticStress(k,l);
		//			end
		//		end
		//		forces(:,ii) = forces(:,ii) + p(:,j) * sum;
        //    
		//	end
		//end
    
		//forces = -volume / 2 * forces;

		double forces[12] = {0}; //rest init to 0.  3X4 matrix.
				
		for (int ii = 0; ii < 4; ii++)
		{
			for (int j = 0; j < 4; j++)
			{
				double sum = 0;
				for (int k = 0; k < 3; k++)
				{
					for (int l = 0; l < 3; l++)
					{
						sum = sum + beta[j * 4 * numTetra + currentTetrad * 4 + l] * beta[ii*4*numTetra + currentTetrad*4 + k] * elasticStress[k * 3 + l];
					}
				}
				for (int k = 0; k < 3; k++)
				{
					forces[k*4+ii] += p[k*4+j] * sum;
				}
            
			}
		}

		//forces = -volume / 2 * forces;
		for (int i = 0; i < 12; i++)
		{
			forces[i] = -volume / 2 * forces[i];
		}
    
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print3By4MatrixSingleIndex(forces,"forces", logger -> MEDIUM);
		}
		#endif


		//currentForce(:, triangles(1, i)) = forces(:,1) - kd * inVelocities(:,triangles(1,i));
		//currentForce(:, triangles(2, i)) = forces(:,2) - kd * inVelocities(:,triangles(2,i));
		//currentForce(:, triangles(3, i)) = forces(:,3) - kd * inVelocities(:,triangles(3,i));
		//currentForce(:, triangles(4, i)) = forces(:,4) - kd * inVelocities(:,triangles(4,i));
    
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				//FIX: = changed to +=
				currentForce[i * numVertices + tetraList[j * numTetra + currentTetrad]] += forces[i*4+j] - kd * defVertices[tetraList[j*numTetra+currentTetrad]].velocity[i];
			}
		}
		
		
	}
	

	//inVelocities = inVelocities + currentForce * elapsedTime + repmat(earthGravity, 1, size(inVelocities,2)) * elapsedTime;
	//defVertices = defVertices + inVelocities * elapsedTime;


	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printVertexTypeMatrix(currentForce,numVertices,"currentForce", logger ->MEDIUM);
	}
	#endif

	

	if (isAnimating)
	{
		for (int i = 0; i < numVertices; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				//Add delta velocity to each particle's velocity
				defVertices[i].velocity[j] += (currentForce[j * numVertices + i] / massMatrix[i]) * deltaT;

				if (j == 1)
				{
					defVertices[i].velocity[j] += -earthGravityValue * deltaT;

				}

				//Use explicit integration with velocity to update each particle's position
				defVertices[i].position[j] += defVertices[i].velocity[j] * deltaT;
			}
		
		}

		doCollisionDetectionAndResponse(deltaT);

		
	}

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger ->printVelocitiesAndPositions(defVertices,numVertices,"Deformed Positions", "Deformed Velocities", logger ->MEDIUM);
	}
	#endif

	timeSinceVideoWrite += deltaT;

	iter++;	

	

	
	
}