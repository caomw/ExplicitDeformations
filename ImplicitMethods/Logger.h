#pragma once

//#include "Edge.h"
//#include "Particle.h"
#include "Vertex.h"
#include <fstream>
#include <string>

using namespace std;

//#define DEBUGGING 1				//If present, allows debugging/logging logic to be compiled.  If absent, it will be omitted, resulting in faster execution.

//This class implements most logging and manages the settings for it
//Currently it only logs to the console.
//Field loggingLevel can be used to increase or decrease the amount of logging
class Logger
{
public:
	//This enum represents differrent levels of logging
	//Each time a logging method is called, the last parameter represents the requested logging level for that method
	//Statements will not log unless the current logging level is at least as high as the requested level
	enum LoggingLevel
	{
		LIGHT = 0,		//Least logging (currently only timing information)
		MEDIUM = 1,		//Normal logging (can still be VERY big on the console unless you use very few particles)
		FULL = 2		//Maximumn logging (logs the conjugate gradient method info)
	};
	bool isLogging;						//true if logging (cout's) should be shown; false if not
	bool printConstraintDeltaV;			//true if deltaV's should be printed for constraint particles (separate from isDebugging setting)

	LoggingLevel loggingLevel;	//Stores the value of the current logging level
	ofstream logFile;


	Logger();
	void initFile(string name);
	//void printAFull(double * ADiagonalEntries, int numParticles, double * ANonDiagonalEntries, int numEdges, Edge * edgeList, char * name, LoggingLevel requestLevel);
	void printVector(double * vector, int size, char * message, LoggingLevel requestLevel);
	void print3By3Matrix(double matrix [3][3], char * message, LoggingLevel requestLevel);
	void print3By3MatrixSingleIndex(double * matrix, int numMatrices, int offSet, char * message, LoggingLevel requestLevel);
	void print3By3MatrixSingleIndex(double matrix [9], char * message, LoggingLevel requestLevel);
	void print3By4Matrix(double matrix [3][4], char * message, LoggingLevel requestLevel);
	void print3By4MatrixSingleIndex(double matrix [12], char * message, LoggingLevel requestLevel);
	void print3By4MatrixSingleIndex(double * matrix, int numMatrices, int offSet, char * message, LoggingLevel requestLevel);
	void print4By4MatrixSingleIndex(double matrix[16], char * message, LoggingLevel requestLevel);
	void print4By4MatrixSingleIndex(double * matrix, int numMatrices, int offSet, char * message, LoggingLevel requestLevel);
	void printNormals(double * normals, int numVertices, char * message, LoggingLevel requestLevel);
	void printStrainToStressMatrix(double strainToStress[6][6], char * message, LoggingLevel requestLevel);
	void printStrainToStressMatrixSingleIndex(double strainToStress[36], char * message, LoggingLevel requestLevel);
	//void printParticles(Particle * particles, int numParticles, LoggingLevel requestLevel);
	//void printNormals(Particle * particles, int numParticles, LoggingLevel requestLevel);
	//void printEdges(Edge * edges, int numEdges, LoggingLevel requestLevel);
	void printMassMatrix(double * massMatrix, int numParticles, LoggingLevel requestLevel);
	void printVertexTypeMatrix(double * matrix, int numVertices, char * message, LoggingLevel requestLevel);
	void printVelocitiesAndPositions(Vertex * vertices, int numVertices, char * message, char * message2, LoggingLevel requestLevel);
	void printText(string text);
	void printIteration(string text, int number);
};
