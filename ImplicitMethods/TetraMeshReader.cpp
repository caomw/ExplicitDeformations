#include <iostream>
#include "TetraMeshReader.h"

using namespace std;

//Meshes from and format based on: http://www.cs.berkeley.edu/~jrs/stellar/#anims

//This method opens the stellar file but does not start reading it
bool TetraMeshReader::openFile(char * nodeFileName, char * elementFileName)
{
	nodeFile.open(nodeFileName, ios::in);
	elementFile.open(elementFileName, ios::in);

	if (!nodeFile.is_open())
	{
		cerr << "Issue loading node file " << nodeFileName << endl;
		return false;
	}
	else
	{
		if (!elementFile.is_open())
		{
			cerr << "Issue loading element file " << elementFileName << endl;
			return false;
		}
		else
		{
			return true;
		}
	}

}

//This method loads the data for a stellar input file
bool TetraMeshReader::loadData(Vertex *& vertexList, int & vertexCount, int *& tetraList, int & tetraCount, Logger * logger)
{
	if (!nodeFile.is_open() || !elementFile.is_open())
	{
		cerr << "Files not opened properly for reading" << endl;
		return false;
	}

	nodeFile >> vertexCount; //# vertices

	if (!nodeFile)
	{
		cerr << "Issue in loading data for node file" << endl;
		return false;
	}

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		if (logger -> loggingLevel >= logger ->MEDIUM)
		{
			cout << "Vertex Count: " << vertexCount << endl;
		}
	}
	#endif

	vertexList = new Vertex[vertexCount];

	int temp;

	for (int i = 1; i <= 3; i++)
	{
		nodeFile >> temp; //dimension, #attributes, #boundary marker
		if (!nodeFile)
		{
			cerr << "Issue in loading data for node file" << endl;
			return false;
		}
	}
	
	int vertexNumber = 0;
	bool hadProblem = true;
	while (nodeFile)
	{
		nodeFile >> temp; //Point number (assuming they're all ordered)

		if (!nodeFile)
		{
			if (nodeFile.eof()) //At the end of the file, this is where it will be caught.  There is no problem for eof.
			{
				hadProblem = false;
			}
			continue;  //Let the normal loop catch file failures and deal with them
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger -> loggingLevel >= logger ->MEDIUM)
			{
				cout << "Vertex #" << temp << ":" << endl;
			}
		}
		#endif

		for (int i = 0; i < 3; i++)
		{
			nodeFile >> vertexList[vertexNumber].position[i];
			if (!nodeFile)
			{
				break;	//Break out of for loop
			}

			vertexList[vertexNumber].velocity[i] = 0;
			//vertexList[vertexNumber].normal[i] = 0;
		}
		if (!nodeFile)
		{
			continue; //Let the loop catch the file failure and eal with it
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger -> loggingLevel >= logger ->MEDIUM)
			{
				//cout << "Read a vertex" << endl;
				for (int i = 0; i < 3; i++)
				{
					cout << vertexList[vertexNumber].position[i] << " ";
				}
				cout << endl;
			}
		}
		#endif

		nodeFile >> temp; //Boundary field

		vertexNumber++;
	}

	if (hadProblem)
	{
		cerr << "Issue in loading data for node file" << endl;
		return false;
	}

	////////////////////////////////////////
	///////////////////////////////////////
	elementFile >> tetraCount; //# tetrahedra

	if (!elementFile)
	{
		cerr << "Issue in loading data for element file" << endl;
		return false;
	}

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		if (logger -> loggingLevel >= logger ->MEDIUM)
		{
			cout << "Tetra Count: " << tetraCount << endl;
		}
	}
	#endif

	tetraList = new int [tetraCount * 4];

	for (int i = 1; i <= 2; i++)
	{
		elementFile >> temp; //points per tetrahedron, #attributes
		if (!elementFile)
		{
			cerr << "Issue in loading data for node file" << endl;
			return false;
		}
	}
	
	int tetraNumber = 0;
	hadProblem = true;
	while (elementFile)
	{
		elementFile >> temp; //Tetrahedron number (assuming they're all ordered)

		if (!elementFile)
		{
			if (elementFile.eof())
			{
				hadProblem = false;
			}
			continue;  //Let the normal loop catch the file failure and deal with it
		}

		for (int i = 0; i < 4; i++)
		{
			elementFile >> tetraList[i * tetraCount + tetraNumber];
			tetraList[i * tetraCount + tetraNumber]--; //Convert from 1 based index system to 0 based index system
			if (!elementFile)
			{
				break;	//Break out of for loop
			}
		}

		if (!elementFile)
		{
			continue;
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger -> loggingLevel >= logger ->MEDIUM)
			{
				cout << "Read tetrahedron:" << endl;
				for (int i = 0; i < 4; i++)
				{
					cout << tetraList[i * tetraCount + tetraNumber] << " ";
				}
				cout << endl;
			}
		}
		#endif

		tetraNumber++;
	} //while (elementFile)

	if (hadProblem)
	{
		cerr << "Issue in loading data for element file" << endl;
		return false;
	}

	return true;


}

//This method closes the stellar file
void TetraMeshReader::closeFile()
{
	nodeFile.close();
	elementFile.close();
}