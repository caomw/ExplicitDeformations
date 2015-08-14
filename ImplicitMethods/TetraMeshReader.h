#pragma once

#include "ParticleSystem.h"
#include <fstream>

using namespace std;

//This class reads a tetrahedral model from a Stellar input file
//Meshes from and format based on: http://www.cs.berkeley.edu/~jrs/stellar/#anims
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

