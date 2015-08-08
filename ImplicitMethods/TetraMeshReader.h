#pragma once

#include "ParticleSystem.h"
#include <fstream>

using namespace std;

//This class reads a tetrahedral model from a Stellar input file
class TetraMeshReader
{
private:
	fstream nodeFile;		
	fstream elementFile;
public:
	bool openFile(char * nodeFileName, char *elementFileName);
	bool loadData(Vertex *& vertexList, int & vertexCount, int *& tetraList, int & tetraCount, Logger * logger);
	void closeFile();
};

