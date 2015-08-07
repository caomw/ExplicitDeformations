#include <iostream>
#include <assert.h>
#include <math.h>
#include <limits>
#include "Logger.h"
#include "StanfordSystem.h"
#include <iomanip>
#include <fstream>

using namespace std;

#include "Macros.h"

const double epsilon = 1e-12;	//Used to check approximate equality to 0
int iteration = 1;

StanfordSystem::StanfordSystem(Vertex * vertexList, int vertexCount, int * tetraList, int tetraCount, Logger * logger) : ParticleSystem(vertexList, vertexCount, tetraList, tetraCount, logger)
{
	strcpy(text, "Method 1");
	double K = 700;
	mu = 700; // change hack
	lambda = K - (2.0/3) * mu;		//Lame's first parameter
	kd = 0.2;

	cout << fixed << setprecision(5);
	logger ->initFile("stanfordDeformation.log");
	logger ->printText(text);

	//Optimization - Precompute cross products and force contribution sums derived from
	//the original vertices
	
	crossProductSums = new double[DIMENSION * numTetra * 4];
	invDm = new double[DIMENSION * DIMENSION * numTetra]; //1 3X3 matrix for each tetrahedron

	//for i = 1:size(triangles,2)
	for (int i = 0; i < numTetra; i++)
	{
		//Precompute the cross products executed on the original mesh (this can be used in place of areas
		//times normals)
		//vectorDifferences = zeros(3,4);
		//vectorDifferences(:, 1) = orgVertices(:, triangles(2, i)) - orgVertices(:, triangles(1,i));
		//vectorDifferences(:, 2) = orgVertices(:, triangles(3, i)) - orgVertices(:, triangles(2,i));
		//vectorDifferences(:, 3) = orgVertices(:, triangles(1, i)) - orgVertices(:, triangles(4,i));
		//vectorDifferences(:, 4) = orgVertices(:, triangles(3, i)) - orgVertices(:, triangles(4,i));

		double vectorDifferences [DIMENSION][4];
		for (int j = 0; j < DIMENSION; j++)
		{
			//FIX - changed defVertices below to orgVertices - was clearly causing problems
			//vectorDifferences[j][0] = defVertices[tetraList[1 * numTetra + i]].position[j] - defVertices[tetraList[0 * numTetra + i]].position[j];
			//vectorDifferences[j][1] = defVertices[tetraList[2 * numTetra + i]].position[j] - defVertices[tetraList[1 * numTetra + i]].position[j];
			//vectorDifferences[j][2] = defVertices[tetraList[0 * numTetra + i]].position[j] - defVertices[tetraList[3 * numTetra + i]].position[j];
			//vectorDifferences[j][3] = defVertices[tetraList[2 * numTetra + i]].position[j] - defVertices[tetraList[3 * numTetra + i]].position[j];
			vectorDifferences[j][0] = orgVertices[tetraList[1 * numTetra + i]].position[j] - orgVertices[tetraList[0 * numTetra + i]].position[j];
			vectorDifferences[j][1] = orgVertices[tetraList[2 * numTetra + i]].position[j] - orgVertices[tetraList[1 * numTetra + i]].position[j];
			vectorDifferences[j][2] = orgVertices[tetraList[0 * numTetra + i]].position[j] - orgVertices[tetraList[3 * numTetra + i]].position[j];
			vectorDifferences[j][3] = orgVertices[tetraList[2 * numTetra + i]].position[j] - orgVertices[tetraList[3 * numTetra + i]].position[j];
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print3By4Matrix(vectorDifferences,"vectorDifferences", logger ->MEDIUM);
		}
		#endif

		//crossProducts(:, (i-1)*4 + 1) = crossProduct(vectorDifferences(:,2), vectorDifferences(:,1)); %Vertices 1 2 and 3;
		//crossProducts(:, (i-1)*4 + 2) = crossProduct(vectorDifferences(:,4), vectorDifferences(:,2)); %Vertices 2 3 and 4;
		//crossProducts(:, (i-1)*4 + 3) = crossProduct(vectorDifferences(:,3), vectorDifferences(:,1)); %Vertices 1 2 and 4;
		//crossProducts(:, (i-1)*4 + 4) = crossProduct(vectorDifferences(:,4), vectorDifferences(:,3)); %Vertices 1 3 and 4;

		//TO DO (TODO): Rename array normals as crossProducts and give it local scope (not class scope) in order to save storage
		normals[0 * 4 * numTetra + (i * 4 + 0)] = (-vectorDifferences[1][1] * vectorDifferences[2][0] + vectorDifferences[2][1] * vectorDifferences[1][0]);	
		normals[1 * 4 * numTetra + (i * 4 + 0)] = (-vectorDifferences[2][1] * vectorDifferences[0][0] + vectorDifferences[0][1] * vectorDifferences[2][0]);	
		normals[2 * 4 * numTetra + (i * 4 + 0)] = (-vectorDifferences[0][1] * vectorDifferences[1][0] + vectorDifferences[1][1] * vectorDifferences[0][0]);	

		normals[0 * 4 * numTetra + (i * 4 + 1)] = (-vectorDifferences[1][3] * vectorDifferences[2][1] + vectorDifferences[2][3] * vectorDifferences[1][1]);	
		normals[1 * 4 * numTetra + (i * 4 + 1)] = (-vectorDifferences[2][3] * vectorDifferences[0][1] + vectorDifferences[0][3] * vectorDifferences[2][1]);	
		normals[2 * 4 * numTetra + (i * 4 + 1)] = (-vectorDifferences[0][3] * vectorDifferences[1][1] + vectorDifferences[1][3] * vectorDifferences[0][1]);	

		normals[0 * 4 * numTetra + (i * 4 + 2)] = (-vectorDifferences[1][2] * vectorDifferences[2][0] + vectorDifferences[2][2] * vectorDifferences[1][0]);	
		normals[1 * 4 * numTetra + (i * 4 + 2)] = (-vectorDifferences[2][2] * vectorDifferences[0][0] + vectorDifferences[0][2] * vectorDifferences[2][0]);	
		normals[2 * 4 * numTetra + (i * 4 + 2)] = (-vectorDifferences[0][2] * vectorDifferences[1][0] + vectorDifferences[1][2] * vectorDifferences[0][0]);	

		normals[0 * 4 * numTetra + (i * 4 + 3)] = (-vectorDifferences[1][3] * vectorDifferences[2][2] + vectorDifferences[2][3] * vectorDifferences[1][2]);	
		normals[1 * 4 * numTetra + (i * 4 + 3)] = (-vectorDifferences[2][3] * vectorDifferences[0][2] + vectorDifferences[0][3] * vectorDifferences[2][2]);	
		normals[2 * 4 * numTetra + (i * 4 + 3)] = (-vectorDifferences[0][3] * vectorDifferences[1][2] + vectorDifferences[1][3] * vectorDifferences[0][2]);	
    
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			//Note: We may not need the old logger's printNormals method anymore
			logger ->print3By4MatrixSingleIndex(normals, numTetra, i, "initial normals", logger ->MEDIUM);
		}
		#endif

		//crossProductSums(:, (i-1)*4 + 1) = -1/6 * (crossProducts(:, (i-1)*4 + 1) + crossProducts(:, (i-1)*4 + 3) + crossProducts(:, (i-1)*4 + 4));
		//crossProductSums(:, (i-1)*4 + 2) = -1/6 * (crossProducts(:, (i-1)*4 + 1) + crossProducts(:, (i-1)*4 + 3) + crossProducts(:, (i-1)*4 + 2));
		//crossProductSums(:, (i-1)*4 + 3) = -1/6 * (crossProducts(:, (i-1)*4 + 1) + crossProducts(:, (i-1)*4 + 2) + crossProducts(:, (i-1)*4 + 4));
		//crossProductSums(:, (i-1)*4 + 4) = -1/6 * (crossProducts(:, (i-1)*4 + 2) + crossProducts(:, (i-1)*4 + 3) + crossProducts(:, (i-1)*4 + 4));

		for (int j = 0; j < DIMENSION; j++)
		{
			crossProductSums[numTetra * 4 * j + i*4 + 0] = -1.0/6 * (normals[numTetra * 4 * j + i*4 + 0] + normals[numTetra * 4 * j + i*4 + 2] + normals[numTetra * 4 * j + i*4 + 3]);
			crossProductSums[numTetra * 4 * j + i*4 + 1] = -1.0/6 * (normals[numTetra * 4 * j + i*4 + 0] + normals[numTetra * 4 * j + i*4 + 2] + normals[numTetra * 4 * j + i*4 + 1]);
			crossProductSums[numTetra * 4 * j + i*4 + 2] = -1.0/6 * (normals[numTetra * 4 * j + i*4 + 0] + normals[numTetra * 4 * j + i*4 + 1] + normals[numTetra * 4 * j + i*4 + 3]);
			crossProductSums[numTetra * 4 * j + i*4 + 3] = -1.0/6 * (normals[numTetra * 4 * j + i*4 + 1] + normals[numTetra * 4 * j + i*4 + 2] + normals[numTetra * 4 * j + i*4 + 3]);
		}

		//%Precompute inverse matrices and other associated matrices
		//%opti - You can store the inverse matrix (and Dm), but you need to use
		//%the memory to store it, for EACH element (tetrahedron)
		//%To do: reevaluate vector subtraction order
    
		//Dm(:, 1) = orgVertices(:, triangles(1, i)) - orgVertices(:, triangles(2,i));
		//Dm(:, 2) = orgVertices(:, triangles(3, i)) - orgVertices(:, triangles(2,i));
		//Dm(:, 3) = orgVertices(:, triangles(4, i)) - orgVertices(:, triangles(2,i));
    
		//if isDebugging == 1
		//	display('Dm is:');  
		//	display(Dm);
		//end
    
		//invDm(:,(i-1)*3 + 1:(i-1)*3 + 3) = inv(Dm);

		//To do: reevaluate vector subtraction order
		//Dm = zeros(3,3);
		double Dm[DIMENSION][DIMENSION];
		for (int j = 0; j < DIMENSION; j++)
		{
			Dm[j][0] = orgVertices[tetraList[0 * numTetra + i]].position[j] - orgVertices[tetraList[1 * numTetra + i]].position[j];
			Dm[j][1] = orgVertices[tetraList[2 * numTetra + i]].position[j] - orgVertices[tetraList[1 * numTetra + i]].position[j];
			Dm[j][2] = orgVertices[tetraList[3 * numTetra + i]].position[j] - orgVertices[tetraList[1 * numTetra + i]].position[j];
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print3By3Matrix(Dm,"Dm", logger ->MEDIUM);
		}
		#endif


		//Find inv(Dm)
		//double invDm[DIMENSION][DIMENSION];
	
		double invDeterminant = 
		Dm[0][0] * (Dm[2][2] * Dm[1][1] - Dm[2][1] * Dm[1][2])
		- Dm[1][0] * (Dm[2][2] * Dm[0][1] - Dm[2][1] * Dm[0][2])
		+ Dm[2][0] * (Dm[1][2] * Dm[0][1] - Dm[1][1] * Dm[0][2]);
		invDeterminant = 1 / invDeterminant;

		//Store the inverse matrix for later use
		invDm[0 * numTetra * DIMENSION + i * DIMENSION + 0] = invDeterminant * (Dm[2][2] * Dm[1][1] - Dm[2][1] * Dm[1][2]);
		invDm[0 * numTetra * DIMENSION + i * DIMENSION + 1] = invDeterminant * (-(Dm[2][2] * Dm[0][1] - Dm[2][1] * Dm[0][2]));
		invDm[0 * numTetra * DIMENSION + i * DIMENSION + 2] = invDeterminant * (Dm[1][2] * Dm[0][1] - Dm[1][1] * Dm[0][2]);
		invDm[1 * numTetra * DIMENSION + i * DIMENSION + 0] = invDeterminant * (-(Dm[2][2] * Dm[1][0] - Dm[2][0] * Dm[1][2]));
		invDm[1 * numTetra * DIMENSION + i * DIMENSION + 1] = invDeterminant * (Dm[2][2] * Dm[0][0] - Dm[2][0] * Dm[0][2]);
		invDm[1 * numTetra * DIMENSION + i * DIMENSION + 2] = invDeterminant * (-(Dm[1][2] * Dm[0][0] - Dm[1][0] * Dm[0][2]));
		invDm[2 * numTetra * DIMENSION + i * DIMENSION + 0] = invDeterminant * (Dm[2][1] * Dm[1][0] - Dm[2][0] * Dm[1][1]);
		invDm[2 * numTetra * DIMENSION + i * DIMENSION + 1] = invDeterminant * (-(Dm[2][1] * Dm[0][0] - Dm[2][0] * Dm[0][1]));
		invDm[2 * numTetra * DIMENSION + i * DIMENSION + 2] = invDeterminant * (Dm[1][1] * Dm[0][0] - Dm[1][0] * Dm[0][1]);
	
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print3By3MatrixSingleIndex(invDm, numTetra, i, "Inverse of DM", logger -> MEDIUM);
		}
		#endif


	//end
	}

	


}

StanfordSystem::~StanfordSystem()
{
	delete [] crossProductSums;
	delete [] invDm;
}

void StanfordSystem::doUpdate(double deltaT)
{
	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printText("=======================================================================================================");
		logger -> printIteration("Iteration #", iteration);
		logger ->printVelocitiesAndPositions(defVertices, numVertices, "Def Positions - start of method", "Velocities - start of method", logger ->MEDIUM);
		logger ->printVelocitiesAndPositions(orgVertices, numVertices, "Org Positions - start of method", "Velocities - start of method", logger ->MEDIUM);
		if (logger -> loggingLevel >= logger ->MEDIUM)
		{
			for (int i = 0; i < numTetra; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					cout << tetraList[j * numTetra + i] << " ";
				}
				cout << endl;
			}
		}
	}
	#endif

	//Start with 0 force each iteration
	for (int i = 0; i < DIMENSION * numVertices; i++)
	{
		currentForce[i] = 0;
	}

	for (int i = 0; i < DIMENSION * numTetra * 4; i++)
	{
		normals[i] = 0;
	}

	for (int i = 0; i < numTetra; i++)
	{
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> printText("**************");
			logger -> printIteration("Tetrahedron #", i);
		}
		#endif
		   
    //display('Dm is:');  
    //display(Dm);

    double Ds [DIMENSION][DIMENSION];
    
	for (int j = 0; j < DIMENSION; j++)
	{
		Ds[j][0] = defVertices[tetraList[0 * numTetra + i]].position[j] - defVertices[tetraList[1 * numTetra + i]].position[j];
		Ds[j][1] = defVertices[tetraList[2 * numTetra + i]].position[j] - defVertices[tetraList[1 * numTetra + i]].position[j];
		Ds[j][2] = defVertices[tetraList[3 * numTetra + i]].position[j] - defVertices[tetraList[1 * numTetra + i]].position[j];
	}

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> print3By3Matrix(Ds,"Ds", logger ->MEDIUM);
	}
	#endif

	//display('Ds is:');
    //display(Ds);

	//F = Ds * inv(Dm);
	double F[DIMENSION][DIMENSION];
	//multiplyMatricesATimesB(F,Ds, DIMENSION, DIMENSION, invDm, DIMENSION, DIMENSION);
	//TO DO (TODO): Look into making a method that multiplies matrices FROM AN OFFSET

	//////////////////////////////////////////////////

	int aRows, bCols, aCols;
	aRows = bCols = aCols = DIMENSION;
	for (int ii = 0; ii < aRows; ii++) /*current row # of calculated final matrix element*/						\
	{																											\
		for (int j = 0; j < bCols; j++) /*current col # of calculated final matrix element*/					\
		{																										\
			F[ii][j] = 0;																						\
			for (int k = 0; k < aCols; k++) /*current col in left matrix and row in right matrix for summing*/	\
			{																									\
				F[ii][j] += Ds[ii][k] * invDm[k * numTetra * DIMENSION + i * DIMENSION + j];														\
			}																									\
																												\
		}																										\
																												\
	}	

	//////////////////////////////////////////////////
	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> print3By3Matrix(F,"F", logger ->MEDIUM);

		//Likely will be used for inverted tetrahedra code:
		//Show us the determinant of F so that we can try to evaluate if it's inverted or not
		//double thatDeterminant = determinant3By3(F);	
	}
	#endif
    
	/*
	double FOneIndex[9];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			FOneIndex[i*3 + j] = F[i][j];
		}
	}

	uninvertF(FOneIndex);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			F[i][j] = FOneIndex[i*3 + j];
		}
	}
	*/

    //greenStrain = (1 / 2) * (F' * F - eye(3));
	//Part I: F' * F
	double greenStrain[DIMENSION][DIMENSION];
	for (int i = 0; i < DIMENSION; i++) //current row # of calculated final matrix element
	{
		for (int j = 0; j < DIMENSION; j++) //current col # of calculated final matrix element
		{
			greenStrain[i][j] = 0;
			for (int k = 0; k < DIMENSION; k++) //current col in left matrix and row in right matrix for summing
			{
				greenStrain[i][j] += F[k][i] * F[k][j]; //We access first F in opposite order since transposed
			}

		}

	}

	//Part 2: Multiply by 0.5 and subtact eye(3) (Identity matrix)
	for (int j = 0; j < DIMENSION; j++)
	{
		for (int k = 0; k < DIMENSION; k++)
		{
			greenStrain[j][k] = 0.5 * (greenStrain[j][k] - (j == k));
		}
	}

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> print3By3Matrix(greenStrain,"Green Strain", logger ->MEDIUM);
	}
	#endif


    //display('The green stress tensor is:');
    //display(greenStrain);

    //Convert stress to force
	//voigtGreenStrain = [greenStrain(1,1) greenStrain(2,2) greenStrain(3,3) ...
    //    2 * greenStrain(2,3) 2 * greenStrain(3,1) 2 * greenStrain(1,2) ]';
	double voigtGreenStrain[6];
	voigtGreenStrain[0] = greenStrain[0][0];
	voigtGreenStrain[1] = greenStrain[1][1];
	voigtGreenStrain[2] = greenStrain[2][2];
	voigtGreenStrain[3] = 2 * greenStrain[1][2];
	voigtGreenStrain[4] = 2 * greenStrain[2][0];
	voigtGreenStrain[5] = 2 * greenStrain[0][1];

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printVector(voigtGreenStrain, 6, "voigtGreenStrain", logger ->MEDIUM);
	}
	#endif
	
	//strainToStress = [ ...
    //    2 * mu + lambda    lambda             lambda              0   0    0;   ...
    //    lambda             2 * mu + lambda    lambda              0   0    0;   ...
    //    lambda             lambda             2 * mu + lambda     0   0    0;   ...
    //    0                   0                 0                   mu  0    0;   ...
    //    0                   0                 0                   0   mu   0;   ...
    //    0                   0                 0                   0   0   mu;   ...
    //];

	double strainToStress[6][6] = {
		{2 * mu + lambda,   lambda         ,   lambda         ,    0,   0,    0},   
		{lambda         ,   2 * mu + lambda,   lambda         ,    0,   0,    0},   
		{lambda         ,   lambda         ,   2 * mu + lambda,    0,   0,    0},   
		{0              ,    0             ,   0              ,    mu,  0,    0},   
		{0              ,    0             ,   0              ,    0,   mu,   0},   
		{0              ,    0             ,   0              ,    0,   0 ,  mu}   
	};

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printStrainToStressMatrix(strainToStress, "strainToStress", logger ->MEDIUM);
	}
	#endif

	
    //voigtStress = strainToStress * voigtGreenStrain;
	double voigtStress[6];
	for (int j = 0; j < 6; j++) //Current row in strainToStress Matrix
	{
		voigtStress[j] = 0;
		for (int k = 0; k < 6; k++) //Current Column in StrainToStressMatrix
		{
			voigtStress[j] += strainToStress[j][k] * voigtGreenStrain[k];
		}
	}

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printVector(voigtStress, 6, "voigtStress", logger ->MEDIUM);
	}
	#endif 

    //Convert our voigt stress vector to the full symmetrix tensor (matrix)
    //secondStress = [
    //    voigtStress(1,1) voigtStress(6,1)    voigtStress(5,1); ...
    //    voigtStress(6,1) voigtStress(2,1)    voigtStress(4,1); ...
    //    voigtStress(5,1) voigtStress(4,1)    voigtStress(3,1)  ...
    //];

	double secondStress[DIMENSION][DIMENSION] = {
		{voigtStress[0], voigtStress[5],    voigtStress[4]},
		{voigtStress[5], voigtStress[1],    voigtStress[3]},
		{voigtStress[4], voigtStress[3],    voigtStress[2]}
	};

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> print3By3Matrix(secondStress,"secondStress", logger ->MEDIUM);
	}
	#endif

	//display('Calclulated the 2nd Piola stress.  It was:');
    //display(secondStress);

    //firstStress = F * secondStress;
	double firstStress[DIMENSION][DIMENSION];
	multiplyMatricesATimesB(firstStress,F,DIMENSION,DIMENSION,secondStress,DIMENSION,DIMENSION);

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> print3By3Matrix(firstStress,"firstStress", logger ->MEDIUM);
	}
	#endif
 
    //currentForce(:,triangles(1,i)) = currentForce(:,triangles(1,i)) + -1/3 * firstStress * (areas(1,1) * normals(:,1) + areas(1,3) * normals(:,3) + areas(1,4) * normals(:,4)) - kd * inVelocities(:,triangles(1,i));
    //currentForce(:,triangles(2,i)) = currentForce(:,triangles(2,i)) + -1/3 * firstStress * (areas(1,1) * normals(:,1) + areas(1,3) * normals(:,3) + areas(1,2) * normals(:,2)) - kd * inVelocities(:,triangles(2,i));
    //currentForce(:,triangles(3,i)) = currentForce(:,triangles(3,i)) + -1/3 * firstStress * (areas(1,1) * normals(:,1) + areas(1,2) * normals(:,2) + areas(1,4) * normals(:,4)) - kd * inVelocities(:,triangles(3,i));
    //currentForce(:,triangles(4,i)) = currentForce(:,triangles(4,i)) + -1/3 * firstStress * (areas(1,2) * normals(:,2) + areas(1,3) * normals(:,3) + areas(1,4) * normals(:,4)) - kd * inVelocities(:,triangles(4,i));


	double temp1[DIMENSION][1];
	double temp2[DIMENSION][1];
	double temp3[DIMENSION][1];
	double temp4[DIMENSION][1];
	for (int j = 0; j < DIMENSION; j++)
	{
		//FIX: added + i * 4 again...
		//temp[j][0] = (areas[0] * normals[j * 4 * numTetra + i * 4 + 0] + areas[2] * normals[j * 4 * numTetra + i * 4 + 2] + areas[3] * normals[j * 4 * numTetra + i * 4 + 3]);
		//temp[j][1] = (areas[0] * normals[j * 4 * numTetra + i * 4 + 0] + areas[2] * normals[j * 4 * numTetra + i * 4 + 2] + areas[1] * normals[j * 4 * numTetra + i * 4 + 1]);
		//temp[j][2] = (areas[0] * normals[j * 4 * numTetra + i * 4 + 0] + areas[1] * normals[j * 4 * numTetra + i * 4 + 1] + areas[3] * normals[j * 4 * numTetra + i * 4 + 3]);
		//temp[j][3] = (areas[1] * normals[j * 4 * numTetra + i * 4 + 1] + areas[2] * normals[j * 4 * numTetra + i * 4 + 2] + areas[3] * normals[j * 4 * numTetra + i * 4 + 3]);
		//TO DO: (TODO) : Get rid of temp - not needed.(?)
		temp1[j][0] = crossProductSums[numTetra*4*j + i*4 + 0];
		temp2[j][0] = crossProductSums[numTetra*4*j + i*4 + 1];
		temp3[j][0] = crossProductSums[numTetra*4*j + i*4 + 2];
		temp4[j][0] = crossProductSums[numTetra*4*j + i*4 + 3];
		
	}

	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger -> printVector(temp1[0],DIMENSION,"temp1", logger ->MEDIUM);
		logger -> printVector(temp2[0],DIMENSION,"temp2", logger ->MEDIUM);
		logger -> printVector(temp3[0],DIMENSION,"temp3", logger ->MEDIUM);
		logger -> printVector(temp4[0],DIMENSION,"temp4", logger ->MEDIUM);
	}
	#endif

	double g1[DIMENSION][1];
	double g2[DIMENSION][1];
	double g3[DIMENSION][1];
	double g4[DIMENSION][1];
	
	multiplyMatricesATimesB(g2,firstStress,DIMENSION,DIMENSION,temp2,DIMENSION, 1);
	multiplyMatricesATimesB(g3,firstStress,DIMENSION,DIMENSION,temp3,DIMENSION, 1);
	multiplyMatricesATimesB(g4,firstStress,DIMENSION,DIMENSION,temp4,DIMENSION, 1);
	for (int i = 0; i < DIMENSION; i++)
	{
		g1[i][0] = -(g2[i][0] + g3[i][0] + g4[i][0]);

	}

	//firstStress * temp   --multiplying the whole firstStress matrix by a matrix with each column being a normal is equivalent to multiplying firstStress matrix by each column separately
		
	for (int j = 0; j < DIMENSION; j++)
	{
		currentForce[j * numVertices + tetraList[0 * numTetra + i]] += g1[j][0] - kd * defVertices[tetraList[0 * numTetra + i]].velocity[j];
		currentForce[j * numVertices + tetraList[1 * numTetra + i]] += g2[j][0] - kd * defVertices[tetraList[1 * numTetra + i]].velocity[j];
		currentForce[j * numVertices + tetraList[2 * numTetra + i]] += g3[j][0] - kd * defVertices[tetraList[2 * numTetra + i]].velocity[j];
		currentForce[j * numVertices + tetraList[3 * numTetra + i]] += g4[j][0] - kd * defVertices[tetraList[3 * numTetra + i]].velocity[j];
	}
    
	//LOGGING
	#ifdef DEBUGGING
	if (logger -> isLogging && logger -> loggingLevel >= logger -> FULL)
	{
		logger ->printText( "Current force for this tetrad in middle of loop");

		logger -> printVertexTypeMatrix(currentForce,numVertices,"currentForce", logger ->MEDIUM);
	}
	#endif
    
    
    
    
	}

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
					defVertices[i].velocity[j] += -earthGravityValue * deltaT; //no mass matrix ref here since earthGravityValue is in fact acceleration
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
	iteration++;

	
}