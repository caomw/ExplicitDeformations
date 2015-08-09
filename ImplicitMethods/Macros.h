//#pragma once
#ifndef MACROS
#define MACROS

#include "Logger.h" //To get DEBUGGING constant

//Macro to find the cross product of two vectors and store the result
//Parameters - resultArray - the name of the array to store the results - must be of size 3 rows X 4 * numberOfTetra elements large (4 normals for each tetrahedron)
//resultIndex - column offset at which to store calcuated normal as a column vector
//arrayA - matrix containing first vector for cross product
//indexA - column at which to find first vector for cross product
//arrayB - matrix containing second vector for cross product
//indexB - column at which to find second vector for cross product
#define crossProduct(resultArray,resultIndex,arrayA, indexA, arrayB, indexB)										\
{																													\
	resultArray[0 * 4 * numTetra + resultIndex] = (arrayA[1][indexA] * arrayB[2][indexB] - arrayA[2][indexA] * arrayB[1][indexB]);	\
	resultArray[1 * 4 * numTetra + resultIndex] = (arrayA[2][indexA] * arrayB[0][indexB] - arrayA[0][indexA] * arrayB[2][indexB]);	\
	resultArray[2 * 4 * numTetra + resultIndex] = (arrayA[0][indexA] * arrayB[1][indexB] - arrayA[1][indexA] * arrayB[0][indexB]);	\
}

//#define crossProductForVertices(resultArray,arrayA, indexA, arrayB, indexB)										\
//{
//    defVertices[tetraList[3 * numTetra + i]].position[1], defVertices[tetraList[3 * numTetra + i]].position[2]
//	resultArray[0 * numVertices + resultIndex] = (arrayA[1][indexA] * arrayB[2][indexB] - arrayA[2][indexA] * arrayB[1][indexB]);	\
//	resultArray[1 * numVertices + resultIndex] = (arrayA[2][indexA] * arrayB[0][indexB] - arrayA[0][indexA] * arrayB[2][indexB]);	\
//	resultArray[2 * numVertices + resultIndex] = (arrayA[0][indexA] * arrayB[1][indexB] - arrayA[1][indexA] * arrayB[0][indexB]);	\
//}


//Macro to multiply two matrices
//Parameters - resultMatrix - matrix in which to store the result (assumed to be size aRows X bCols)
//matrixA - first matrix being multiplied
//aRows - rows in matrixA
//aCols - cols in matrixA
//matrixB - second matrix being multiplied
//bRows - rows in matrixB
//bCols - cols in matrixB
#define multiplyMatricesATimesB(resultMatrix, matrixA, aRows, aCols, matrixB, bRows, bCols)						\
{																												\
	/*Temporary*/																								\
	if (aCols != bRows)																							\
	{																											\
		cerr << "Matrix dimensions don't match -- illegal matrix multiplication." << endl;						\
		assert (aCols == bRows);																				\
	}																											\
																												\
	for (int i = 0; i < aRows; i++) /*current row # of calculated final matrix element*/						\
	{																											\
		for (int j = 0; j < bCols; j++) /*current col # of calculated final matrix element*/					\
		{																										\
			resultMatrix[i][j] = 0;																				\
			for (int k = 0; k < aCols; k++) /*current col in left matrix and row in right matrix for summing*/	\
			{																									\
				resultMatrix[i][j] += matrixA[i][k] * matrixB[k][j];											\
			}																									\
																												\
		}																										\
																												\
	}																											\
																												\
}

//Macro to multiply two matrices - without doubly suscripted arrays
//Parameters - resultMatrix - matrix in which to store the result (assumed to be size aRows X bCols)
//matrixA - first matrix being multiplied
//aRows - rows in matrixA
//aCols - cols in matrixA
//matrixB - second matrix being multiplied
//bRows - rows in matrixB
//bCols - cols in matrixB
//Multiplying (aRows X aCols) * (bRows X bCols) gives us size (aRows X bCols)
#define multiplyMatricesATimesBSingleIndexes(resultMatrix, matrixA, aRows, aCols, matrixB, bRows, bCols)						\
{																												\
	/*Temporary*/																								\
	if (aCols != bRows)																							\
	{																											\
		cerr << "Matrix dimensions don't match -- illegal matrix multiplication." << endl;						\
		assert (aCols == bRows);																				\
	}																											\
																												\
	for (int i = 0; i < aRows; i++) /*current row # of calculated final matrix element*/						\
	{																											\
		for (int j = 0; j < bCols; j++) /*current col # of calculated final matrix element*/					\
		{																										\
			resultMatrix[i * bCols + j] = 0;																				\
			for (int k = 0; k < aCols; k++) /*current col in left matrix and row in right matrix for summing*/	\
			{																									\
				resultMatrix[i * bCols + j] += matrixA[i * aCols + k] * matrixB[k * bCols + j];											\
			}																									\
																												\
		}																										\
																												\
	}																											\
																												\
}

#define crossProductGeneral(resultArray,arrayA, arrayB)														\
{																																	\
	resultArray[0] = (arrayA[1] * arrayB[2] - arrayA[2] * arrayB[1]);	\
	resultArray[1] = (arrayA[2] * arrayB[0] - arrayA[0] * arrayB[2]);	\
	resultArray[2] = (arrayA[0] * arrayB[1] - arrayA[1] * arrayB[0]);	\
}

#define dot3(vector1, vector2) (vector1[0] * vector2[0] + vector1[1] * vector2[1] + vector1[2] * vector2[2])

#define trace3By3(theMatrix) theMatrix[0] + theMatrix[4] + theMatrix[8]

//#define determinant3By3(inMatrix) \
	( \
	  inMatrix[0] * inMatrix[4] * inMatrix[8] + \
	  inMatrix[1] * inMatrix[5] * inMatrix[6] + \
	  inMatrix[2] * inMatrix[3] * inMatrix[7] - \
	  inMatrix[2] * inMatrix[4] * inMatrix[6] - \
	  inMatrix[0] * inMatrix[5] * inMatrix[7] - \
	  inMatrix[1] * inMatrix[3] * inMatrix[8] \
	)

#define determinant3By3(inMatrix) \
	( \
	  inMatrix[0][0] * inMatrix[1][1] * inMatrix[2][2] + \
	  inMatrix[0][1] * inMatrix[1][2] * inMatrix[2][0] + \
	  inMatrix[0][2] * inMatrix[1][0] * inMatrix[2][1] - \
	  inMatrix[0][2] * inMatrix[1][1] * inMatrix[2][0] - \
	  inMatrix[0][0] * inMatrix[1][2] * inMatrix[2][1] - \
	  inMatrix[0][1] * inMatrix[1][0] * inMatrix[2][2] \
	)

#define determinant3By3OneBasedIndices(inMatrix) \
	( \
	  inMatrix[1][1] * inMatrix[2][2] * inMatrix[3][3] + \
	  inMatrix[1][2] * inMatrix[2][3] * inMatrix[3][1] + \
	  inMatrix[1][3] * inMatrix[2][1] * inMatrix[3][2] - \
	  inMatrix[1][3] * inMatrix[2][2] * inMatrix[3][1] - \
	  inMatrix[1][1] * inMatrix[2][3] * inMatrix[3][2] - \
	  inMatrix[1][2] * inMatrix[2][1] * inMatrix[3][3] \
	)

#define determinant3By3SingleIndex(inMatrix) \
	( \
	  inMatrix[0*3+0] * inMatrix[1*3+1] * inMatrix[2*3+2] + \
	  inMatrix[0*3+1] * inMatrix[1*3+2] * inMatrix[2*3+0] + \
	  inMatrix[0*3+2] * inMatrix[1*3+0] * inMatrix[2*3+1] - \
	  inMatrix[0*3+2] * inMatrix[1*3+1] * inMatrix[2*3+0] - \
	  inMatrix[0*3+0] * inMatrix[1*3+2] * inMatrix[2*3+1] - \
	  inMatrix[0*3+1] * inMatrix[1*3+0] * inMatrix[2*3+2] \
	)

#define determinant4By4(inMatrix) \
	( \
	  inMatrix[0*4+0]*inMatrix[1*4+1]*inMatrix[2*4+2]*inMatrix[3*4+3] + inMatrix[0*4+0]*inMatrix[1*4+2]*inMatrix[2*4+3]*inMatrix[3*4+1] + inMatrix[0*4+0]*inMatrix[1*4+3]*inMatrix[2*4+1]*inMatrix[3*4+2] \
	+ inMatrix[0*4+1]*inMatrix[1*4+0]*inMatrix[2*4+3]*inMatrix[3*4+2] + inMatrix[0*4+1]*inMatrix[1*4+2]*inMatrix[2*4+0]*inMatrix[3*4+3] + inMatrix[0*4+1]*inMatrix[1*4+3]*inMatrix[2*4+2]*inMatrix[3*4+0] \
	+ inMatrix[0*4+2]*inMatrix[1*4+0]*inMatrix[2*4+1]*inMatrix[3*4+3] + inMatrix[0*4+2]*inMatrix[1*4+1]*inMatrix[2*4+3]*inMatrix[3*4+0] + inMatrix[0*4+2]*inMatrix[1*4+3]*inMatrix[2*4+0]*inMatrix[3*4+1] \
	+ inMatrix[0*4+3]*inMatrix[1*4+0]*inMatrix[2*4+2]*inMatrix[3*4+1] + inMatrix[0*4+3]*inMatrix[1*4+1]*inMatrix[2*4+0]*inMatrix[3*4+2] + inMatrix[0*4+3]*inMatrix[1*4+2]*inMatrix[2*4+1]*inMatrix[3*4+0] \
	- inMatrix[0*4+0]*inMatrix[1*4+1]*inMatrix[2*4+3]*inMatrix[3*4+2] - inMatrix[0*4+0]*inMatrix[1*4+2]*inMatrix[2*4+1]*inMatrix[3*4+3] - inMatrix[0*4+0]*inMatrix[1*4+3]*inMatrix[2*4+2]*inMatrix[3*4+1] \
	- inMatrix[0*4+1]*inMatrix[1*4+0]*inMatrix[2*4+2]*inMatrix[3*4+3] - inMatrix[0*4+1]*inMatrix[1*4+2]*inMatrix[2*4+3]*inMatrix[3*4+0] - inMatrix[0*4+1]*inMatrix[1*4+3]*inMatrix[2*4+0]*inMatrix[3*4+2] \
	- inMatrix[0*4+2]*inMatrix[1*4+0]*inMatrix[2*4+3]*inMatrix[3*4+1] - inMatrix[0*4+2]*inMatrix[1*4+1]*inMatrix[2*4+0]*inMatrix[3*4+3] - inMatrix[0*4+2]*inMatrix[1*4+3]*inMatrix[2*4+1]*inMatrix[3*4+0] \
	- inMatrix[0*4+3]*inMatrix[1*4+0]*inMatrix[2*4+1]*inMatrix[3*4+2] - inMatrix[0*4+3]*inMatrix[1*4+1]*inMatrix[2*4+2]*inMatrix[3*4+0] - inMatrix[0*4+3]*inMatrix[1*4+2]*inMatrix[2*4+0]*inMatrix[3*4+1] \
	)

void find4By4MatrixInverse(double inMatrix[16], double resultMatrix[16], Logger * logger);

#endif