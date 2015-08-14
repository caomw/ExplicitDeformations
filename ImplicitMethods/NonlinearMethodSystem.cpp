#include <iostream>
#include <iomanip>
#include <assert.h>
#include "Logger.h"
#include "NonlinearMethodSystem.h"
#include "Macros.h"

using namespace std;

const double epsilon = 1e-12;	//Used to check approximate equality to 0

//Based on the paper at http://www-ljk.imag.fr/Publications/Basilic/com.lmc.publi.PUBLI_Article@11f6a0378d9_18c74/tensile.pdf – Simple, yet Accurate Nonlinear Tensile Stiffness
//Pascal Volino et. al

//Constructor
NonlinearMethodSystem::NonlinearMethodSystem(Vertex * vertexList, int vertexCount, int * tetraList, int tetraCount, Logger * logger) : ParticleSystem(vertexList, vertexCount, tetraList, tetraCount, logger)
{
	strcpy(text, "Method 3");
	double K = 600;					//Bulk Modulus
	mu = 600;						//Shear modulus (Lame's second parameter)
	lambda = K - (2.0/3) * mu;		//Lame's first parameter
	kd = 0.25;

	logger -> initFile("nonLinearDeformation.log");
	logger ->printText(text);

	ruWeights = new double [4 * numTetra];
	rvWeights = new double [4 * numTetra];
	rwWeights = new double [4 * numTetra];

	for (int currentTetrad = 0; currentTetrad < numTetra; currentTetrad++)
	{
		//vertex1 = orgVertices(:,triangles(1,currentTriangle));
        //vertex2 = orgVertices(:,triangles(2,currentTriangle));
        //vertex3 = orgVertices(:,triangles(3,currentTriangle));
        //vertex4 = orgVertices(:,triangles(4,currentTriangle));

		//ua = vertex1(1,1);
        //ub = vertex2(1,1);
        //uc = vertex3(1,1);
        //ud = vertex4(1,1);
		double ua = orgVertices[tetraList[0 * numTetra + currentTetrad]].position[0];
		double ub = orgVertices[tetraList[1 * numTetra + currentTetrad]].position[0];
		double uc = orgVertices[tetraList[2 * numTetra + currentTetrad]].position[0];
		double ud = orgVertices[tetraList[3 * numTetra + currentTetrad]].position[0];
        
        //va = vertex1(2,1);
        //vb = vertex2(2,1);
        //vc = vertex3(2,1);
        //vd = vertex4(2,1);
		double va = orgVertices[tetraList[0 * numTetra + currentTetrad]].position[1];
		double vb = orgVertices[tetraList[1 * numTetra + currentTetrad]].position[1];
		double vc = orgVertices[tetraList[2 * numTetra + currentTetrad]].position[1];
		double vd = orgVertices[tetraList[3 * numTetra + currentTetrad]].position[1];
        
        //wa = vertex1(3,1);
        //wb = vertex2(3,1);
        //wc = vertex3(3,1);
        //wd = vertex4(3,1);
		double wa = orgVertices[tetraList[0 * numTetra + currentTetrad]].position[2];
		double wb = orgVertices[tetraList[1 * numTetra + currentTetrad]].position[2];
		double wc = orgVertices[tetraList[2 * numTetra + currentTetrad]].position[2];
		double wd = orgVertices[tetraList[3 * numTetra + currentTetrad]].position[2];

		//ruWeights = [ua ub uc ud; va vb vc vd; wa wb wc wd; 1 1 1 1] \ [1 0 0 0]';
        //rvWeights = [ua ub uc ud; va vb vc vd; wa wb wc wd; 1 1 1 1] \ [0 1 0 0]';
        //rwWeights = [ua ub uc ud; va vb vc vd; wa wb wc wd; 1 1 1 1] \ [0 0 1 0]';
        
		////////////////////////////////////////////////////////////////////////////
		//We need to solve each of the 3 systems above.
		//Since we already have working 4x4 matrix inverse code, it makes perfect sense to use it to solve this system
		//ie:
		//Solve [rua rub ruc rud]' = [ua ub uc ud; va vb vc vd; wa wb wc wd; 1 1 1 1] \ [1 0 0 0]';
		//[rua rub ruc rud]' = inv([ua ub uc ud; va vb vc vd; wa wb wc wd; 1 1 1 1]) * [1 0 0 0]'
		/////////////////////////////////////////////////////////////////////////////////////////

		double squareMatrix[16];
		double inverseSquareMatrix[16];
		
		squareMatrix[0] = ua;   squareMatrix[1] = ub;   squareMatrix[2] = uc;   squareMatrix[3] = ud;
		squareMatrix[4] = va;   squareMatrix[5] = vb;   squareMatrix[6] = vc;   squareMatrix[7] = vd;
		squareMatrix[8] = wa;   squareMatrix[9] = wb;   squareMatrix[10] = wc;  squareMatrix[11] = wd;
		squareMatrix[12] = 1;   squareMatrix[13] = 1;   squareMatrix[14] = 1;   squareMatrix[15] = 1;
		
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print4By4MatrixSingleIndex(squareMatrix,"original vertices square matrix", logger ->MEDIUM);
		}
		#endif

		find4By4MatrixInverse(squareMatrix,inverseSquareMatrix, logger);

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print4By4MatrixSingleIndex(inverseSquareMatrix,"Inverse of original vertices square matrix", logger ->MEDIUM);
		}
		#endif

		double ruWeightsTemp[4];
		double ruRightVector[4] = {1, 0, 0, 0};  //TODO: Can we optimize this matrix-vector multiplication to exploit the zeroes?
		multiplyMatricesATimesBSingleIndexes(ruWeightsTemp,inverseSquareMatrix, 4, 4, ruRightVector, 4, 1);
		
		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger ->printVector(ruWeightsTemp,4,"ru weights", logger ->MEDIUM);
		}
		#endif

		double rvWeightsTemp[4];
		double rvRightVector[4] = {0, 1, 0, 0};  //TODO: Can we optimize this matrix-vector multiplication to exploit the zeroes?
		multiplyMatricesATimesBSingleIndexes(rvWeightsTemp,inverseSquareMatrix, 4, 4, rvRightVector, 4, 1);

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger ->printVector(rvWeightsTemp,4,"rv weights", logger ->MEDIUM);
		}
		#endif

		double rwWeightsTemp[4];
		double rwRightVector[4] = {0, 0, 1, 0};  //TODO: Can we optimize this matrix-vector multiplication to exploit the zeroes?
		multiplyMatricesATimesBSingleIndexes(rwWeightsTemp,inverseSquareMatrix, 4, 4, rwRightVector, 4, 1);

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger ->printVector(rwWeightsTemp,4,"rw weights", logger ->MEDIUM);
		}
		#endif

		for (int i = 0; i < 4; i++)
		{
			ruWeights[currentTetrad*4 + i] = ruWeightsTemp[i];
			rvWeights[currentTetrad*4 + i] = rvWeightsTemp[i];
			rwWeights[currentTetrad*4 + i] = rwWeightsTemp[i];
		}


	}
		
}

NonlinearMethodSystem::~NonlinearMethodSystem()
{
	delete [] ruWeights;
	delete [] rvWeights;
	delete [] rwWeights;
}

//Overridden update method
void NonlinearMethodSystem::doUpdate(double deltaT)
{
	#ifdef DEBUGGING
	if (logger -> isLogging)
	{
		logger ->printVelocitiesAndPositions(defVertices, numVertices, "Positions - start of method", "Velocities - start of method", logger ->MEDIUM);
	}
	#endif
	
    //currentForce = zeros(size(inVelocities,1), size(inVelocities,2));
    
	//Start with 0 force each iteration
	for (int i = 0; i < DIMENSION * numVertices; i++)
	{
		currentForce[i] = 0;
	}

	//for currentTriangle = 1:size(triangles,2)
	for (int currentTetrad = 0; currentTetrad < numTetra; currentTetrad++)
	{
		//rua = ruWeights(1,1);
        //rub = ruWeights(2,1);
        //ruc = ruWeights(3,1);
        //rud = ruWeights(4,1);
        
        //rva = rvWeights(1,1);
        //rvb = rvWeights(2,1);
        //rvc = rvWeights(3,1);
        //rvd = rvWeights(4,1);
                
        //rwa = rwWeights(1,1);
        //rwb = rwWeights(2,1);
        //rwc = rwWeights(3,1);
        //rwd = rwWeights(4,1);
                   
        //vertex1 = defVertices(:,triangles(1,currentTriangle));
        //vertex2 = defVertices(:,triangles(2,currentTriangle));
        //vertex3 = defVertices(:,triangles(3,currentTriangle));
        //vertex4 = defVertices(:,triangles(4,currentTriangle));
        
		//U = rua * vertex1 + rub * vertex2 + ruc * vertex3 + rud * vertex4;
        //V = rva * vertex1 + rvb * vertex2 + rvc * vertex3 + rvd * vertex4;
        //W = rwa * vertex1 + rwb * vertex2 + rwc * vertex3 + rwd * vertex4;

		double U[3], V[3], W[3];
		for (int i = 0; i < DIMENSION; i++)
		{
			U[i] = ruWeights[currentTetrad * 4 + 0] * defVertices[tetraList[0 * numTetra + currentTetrad]].position[i] + ruWeights[currentTetrad * 4 + 1] * defVertices[tetraList[1 * numTetra + currentTetrad]].position[i] + ruWeights[currentTetrad * 4 + 2] * defVertices[tetraList[2 * numTetra + currentTetrad]].position[i] + ruWeights[currentTetrad * 4 + 3] * defVertices[tetraList[3 * numTetra + currentTetrad]].position[i];
			V[i] = rvWeights[currentTetrad * 4 + 0] * defVertices[tetraList[0 * numTetra + currentTetrad]].position[i] + rvWeights[currentTetrad * 4 + 1] * defVertices[tetraList[1 * numTetra + currentTetrad]].position[i] + rvWeights[currentTetrad * 4 + 2] * defVertices[tetraList[2 * numTetra + currentTetrad]].position[i] + rvWeights[currentTetrad * 4 + 3] * defVertices[tetraList[3 * numTetra + currentTetrad]].position[i];
			W[i] = rwWeights[currentTetrad * 4 + 0] * defVertices[tetraList[0 * numTetra + currentTetrad]].position[i] + rwWeights[currentTetrad * 4 + 1] * defVertices[tetraList[1 * numTetra + currentTetrad]].position[i] + rwWeights[currentTetrad * 4 + 2] * defVertices[tetraList[2 * numTetra + currentTetrad]].position[i] + rwWeights[currentTetrad * 4 + 3] * defVertices[tetraList[3 * numTetra + currentTetrad]].position[i];
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger ->printVector(U,3,"U Axis Vector", logger ->MEDIUM);
			logger ->printVector(V,3,"V Axis Vector", logger ->MEDIUM);
			logger ->printVector(W,3,"W Axis Vector", logger ->MEDIUM);
		}
		#endif
        
		//Euu = 0.5 * (U' * U - 1);
        //Evv = 0.5 * (V' * V - 1);
        //Eww = 0.5 * (W' * W - 1);
        //Evw = 0.5 * (V' * W);
        //Euw = 0.5 * (U' * W);
        //Euv = 0.5 * (U' * V);

		double Euu = 0.5 * (dot3(U, U) - 1);
        double Evv = 0.5 * (dot3(V, V) - 1);
        double Eww = 0.5 * (dot3(W, W) - 1);
        double Evw = 0.5 * (dot3(V, W));
        double Euw = 0.5 * (dot3(U, W));
        double Euv = 0.5 * (dot3(U, V));

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger ->loggingLevel >= logger -> MEDIUM)
			{
				cout << "Euu: " << Euu << endl;
				cout << "Evv: " << Evv << endl;
				cout << "Eww: " << Eww << endl;
				cout << "Evw: " << Evw << endl;
				cout << "Euw: " << Euw << endl;
				cout << "Euv: " << Euv << endl;
			}
		}
		#endif

		//3D Hooke's Law based on: https://en.wikipedia.org/w/index.php?title=Hooke%27s_law&oldid=673340434

		//voigtStrain = [Euu Evv Eww 2 * Evw 2 * Euw 2 * Euv]';

		double voigtGreenStrain[6] = {Euu, Evv, Eww, 2 * Evw, 2 * Euw, 2 * Euv};    

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger ->printVector(voigtGreenStrain,6,"voigtGreenStrain", logger ->MEDIUM);
		}
		#endif

        
        //strainToStress = [ ...
        //2 * mu + lambda    lambda             lambda              0   0    0;   ...
        //lambda             2 * mu + lambda    lambda              0   0    0;   ...
        //lambda             lambda             2 * mu + lambda     0   0    0;   ...
        //0                   0                 0                   mu  0    0;   ...
        //0                   0                 0                   0   mu   0;   ...
        //0                   0                 0                   0   0   mu;   ...
        //];

		double strainToStress[36] = 
		{
			2 * mu + lambda,     lambda,             lambda,              0,   0,    0,
			lambda,              2 * mu + lambda,    lambda,              0,   0,    0,
			lambda,              lambda,             2 * mu + lambda,     0,   0,    0,
			0,                   0,                  0,                   mu,  0,    0,
			0,                   0,                  0,                   0,   mu,   0,
			0,                   0,                  0,                   0,   0,    mu
		};

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger ->printStrainToStressMatrixSingleIndex(strainToStress,"strain to stress matrix", logger ->MEDIUM);
		}
		#endif
    
        //voigtStress = strainToStress * voigtStrain;
		double voigtStress[6];

		multiplyMatricesATimesBSingleIndexes(voigtStress, strainToStress,6, 6, voigtGreenStrain, 6,1);

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger ->printVector(voigtStress,6,"voigtStress", logger ->MEDIUM);
		}
		#endif
   
        //ouu = voigtStress(1,1);
        //ovv = voigtStress(2,1);
        //oww = voigtStress(3,1);
        //ovw = voigtStress(4,1);
        //ouw = voigtStress(5,1);
        //ouv = voigtStress(6,1);

		double ouu = voigtStress[0];
        double ovv = voigtStress[1];
        double oww = voigtStress[2];
        double ovw = voigtStress[3];
        double ouw = voigtStress[4];
        double ouv = voigtStress[5];

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger ->loggingLevel >= logger -> MEDIUM)
			{
				logger -> printIteration("ouu:", ouu);
				logger -> printIteration("ovv:", ovv);
				logger -> printIteration("oww:", oww);
				logger -> printIteration("ovw:", ovw);
				logger -> printIteration("ouw:", ouw);
				logger -> printIteration("ouv:", ouv);
			}
		}
		#endif
        
        //volumeMatrix = [ vertex1(1,1) vertex2(1,1) vertex3(1,1) vertex4(1,1); ...
        //                 vertex1(2,1) vertex2(2,1) vertex3(2,1) vertex4(2,1); ...
        //                 vertex1(3,1) vertex2(3,1) vertex3(3,1) vertex4(3,1); ... 
        //                 1            1            1            1             ...
        //               ];

		double volumeMatrix[16] = 
		{ 
			                    defVertices[tetraList[0 * numTetra + currentTetrad]].position[0], defVertices[tetraList[1 * numTetra + currentTetrad]].position[0], defVertices[tetraList[2 * numTetra + currentTetrad]].position[0], defVertices[tetraList[3 * numTetra + currentTetrad]].position[0],
								defVertices[tetraList[0 * numTetra + currentTetrad]].position[1], defVertices[tetraList[1 * numTetra + currentTetrad]].position[1], defVertices[tetraList[2 * numTetra + currentTetrad]].position[1], defVertices[tetraList[3 * numTetra + currentTetrad]].position[1],
								defVertices[tetraList[0 * numTetra + currentTetrad]].position[2], defVertices[tetraList[1 * numTetra + currentTetrad]].position[2], defVertices[tetraList[2 * numTetra + currentTetrad]].position[2], defVertices[tetraList[3 * numTetra + currentTetrad]].position[2],
								                                                               1,                                                                1,                                                                1,                                                                1
		};

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> print4By4MatrixSingleIndex(volumeMatrix,"volume matrix", logger ->MEDIUM);
		}
		#endif
        
		//volume = (1/6) * abs(det(volumeMatrix));
        double volume = (1.0/6) * abs(determinant4By4(volumeMatrix));

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			if (logger ->loggingLevel >= logger -> MEDIUM)
			{
				logger -> printIteration("volume:", volume);
			}
		}
		#endif
		           
        //currentForce(:,triangles(1,currentTriangle)) = currentForce(:,triangles(1,currentTriangle)) - volume * (ouu*rua * U + ouv*(0.5*rua*V + 0.5*rva * U) + ouw*(0.5 * rua*W+0.5*rwa*U) + ovv*rva*V + ovw*(0.5*rva*W + 0.5*rwa*V) + oww*rwa*W) - kd * inVelocities(:,triangles(1,currentTriangle));
        //currentForce(:,triangles(2,currentTriangle)) = currentForce(:,triangles(2,currentTriangle)) - volume * (ouu*rub * U + ouv*(0.5*rub*V + 0.5*rvb * U) + ouw*(0.5 * rub*W+0.5*rwb*U) + ovv*rvb*V + ovw*(0.5*rvb*W + 0.5*rwb*V) + oww*rwb*W) - kd * inVelocities(:,triangles(2,currentTriangle));
        //currentForce(:,triangles(3,currentTriangle)) = currentForce(:,triangles(3,currentTriangle)) - volume * (ouu*ruc * U + ouv*(0.5*ruc*V + 0.5*rvc * U) + ouw*(0.5 * ruc*W+0.5*rwc*U) + ovv*rvc*V + ovw*(0.5*rvc*W + 0.5*rwc*V) + oww*rwc*W) - kd * inVelocities(:,triangles(3,currentTriangle));
        //currentForce(:,triangles(4,currentTriangle)) = currentForce(:,triangles(4,currentTriangle)) - volume * (ouu*rud * U + ouv*(0.5*rud*V + 0.5*rvd * U) + ouw*(0.5 * rud*W+0.5*rwd*U) + ovv*rvd*V + ovw*(0.5*rvd*W + 0.5*rwd*V) + oww*rwd*W) - kd * inVelocities(:,triangles(4,currentTriangle));
				
		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < DIMENSION; i++)
			{
				currentForce[i * numVertices + tetraList[j * numTetra + currentTetrad]] +=  -volume * (ouu*ruWeights[currentTetrad * 4 + j] * U[i] + ouv*(0.5*ruWeights[currentTetrad * 4 + j]*V[i] + 0.5*rvWeights[currentTetrad * 4 + j] * U[i]) + ouw*(0.5 * ruWeights[currentTetrad * 4 + j]*W[i]+0.5*rwWeights[currentTetrad * 4 + j]*U[i]) + ovv*rvWeights[currentTetrad * 4 + j]*V[i] + ovw*(0.5*rvWeights[currentTetrad * 4 + j]*W[i] + 0.5*rwWeights[currentTetrad * 4 + j]*V[i]) + oww*rwWeights[currentTetrad * 4 + j]*W[i]) - kd * defVertices[tetraList[j*numTetra+currentTetrad]].velocity[i];
			}
		}

		#ifdef DEBUGGING
		if (logger -> isLogging)
		{
			logger -> printVertexTypeMatrix(currentForce,numVertices,"currentForce", logger ->MEDIUM);
		}
		#endif
        
    //end
	}
	    
    //currentForce = currentForce + repmat(earthGravity, 1, size(inVelocities,2));
    //inVelocities = inVelocities + currentForce * elapsedTime;
    
    //defVertices = defVertices + inVelocities * elapsedTime;

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
    
    timeSinceVideoWrite += deltaT;
}

