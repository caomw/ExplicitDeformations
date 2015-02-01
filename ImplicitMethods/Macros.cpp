#include "Macros.h"
#include "Logger.h"

//This method finds the inverse of a 4 by 4 matrix (inMatrix)
//It is output to resultMatrix
//Logger is a reference to the logger (not currently used)
void find4By4MatrixInverse(double inMatrix[16], double resultMatrix[16], Logger * logger)
{
	//From http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
	double invDeterminant =
    1.0 / 
	(
	  inMatrix[0*4+0]*inMatrix[1*4+1]*inMatrix[2*4+2]*inMatrix[3*4+3] + inMatrix[0*4+0]*inMatrix[1*4+2]*inMatrix[2*4+3]*inMatrix[3*4+1] + inMatrix[0*4+0]*inMatrix[1*4+3]*inMatrix[2*4+1]*inMatrix[3*4+2]
	+ inMatrix[0*4+1]*inMatrix[1*4+0]*inMatrix[2*4+3]*inMatrix[3*4+2] + inMatrix[0*4+1]*inMatrix[1*4+2]*inMatrix[2*4+0]*inMatrix[3*4+3] + inMatrix[0*4+1]*inMatrix[1*4+3]*inMatrix[2*4+2]*inMatrix[3*4+0]
	+ inMatrix[0*4+2]*inMatrix[1*4+0]*inMatrix[2*4+1]*inMatrix[3*4+3] + inMatrix[0*4+2]*inMatrix[1*4+1]*inMatrix[2*4+3]*inMatrix[3*4+0] + inMatrix[0*4+2]*inMatrix[1*4+3]*inMatrix[2*4+0]*inMatrix[3*4+1]
	+ inMatrix[0*4+3]*inMatrix[1*4+0]*inMatrix[2*4+2]*inMatrix[3*4+1] + inMatrix[0*4+3]*inMatrix[1*4+1]*inMatrix[2*4+0]*inMatrix[3*4+2] + inMatrix[0*4+3]*inMatrix[1*4+2]*inMatrix[2*4+1]*inMatrix[3*4+0]
	- inMatrix[0*4+0]*inMatrix[1*4+1]*inMatrix[2*4+3]*inMatrix[3*4+2] - inMatrix[0*4+0]*inMatrix[1*4+2]*inMatrix[2*4+1]*inMatrix[3*4+3] - inMatrix[0*4+0]*inMatrix[1*4+3]*inMatrix[2*4+2]*inMatrix[3*4+1]
	- inMatrix[0*4+1]*inMatrix[1*4+0]*inMatrix[2*4+2]*inMatrix[3*4+3] - inMatrix[0*4+1]*inMatrix[1*4+2]*inMatrix[2*4+3]*inMatrix[3*4+0] - inMatrix[0*4+1]*inMatrix[1*4+3]*inMatrix[2*4+0]*inMatrix[3*4+2]
	- inMatrix[0*4+2]*inMatrix[1*4+0]*inMatrix[2*4+3]*inMatrix[3*4+1] - inMatrix[0*4+2]*inMatrix[1*4+1]*inMatrix[2*4+0]*inMatrix[3*4+3] - inMatrix[0*4+2]*inMatrix[1*4+3]*inMatrix[2*4+1]*inMatrix[3*4+0]
	- inMatrix[0*4+3]*inMatrix[1*4+0]*inMatrix[2*4+1]*inMatrix[3*4+2] - inMatrix[0*4+3]*inMatrix[1*4+1]*inMatrix[2*4+2]*inMatrix[3*4+0] - inMatrix[0*4+3]*inMatrix[1*4+2]*inMatrix[2*4+0]*inMatrix[3*4+1]
	);

	resultMatrix[0*4+0] = invDeterminant * (inMatrix[1*4+1]*inMatrix[2*4+2]*inMatrix[3*4+3] + inMatrix[1*4+2]*inMatrix[2*4+3]*inMatrix[3*4+1] + inMatrix[1*4+3]*inMatrix[2*4+1]*inMatrix[3*4+2] - inMatrix[1*4+1]*inMatrix[2*4+3]*inMatrix[3*4+2] - inMatrix[1*4+2]*inMatrix[2*4+1]*inMatrix[3*4+3] - inMatrix[1*4+3]*inMatrix[2*4+2]*inMatrix[3*4+1]);
	resultMatrix[0*4+1] = invDeterminant * (inMatrix[0*4+1]*inMatrix[2*4+3]*inMatrix[3*4+2] + inMatrix[0*4+2]*inMatrix[2*4+1]*inMatrix[3*4+3] + inMatrix[0*4+3]*inMatrix[2*4+2]*inMatrix[3*4+1] - inMatrix[0*4+1]*inMatrix[2*4+2]*inMatrix[3*4+3] - inMatrix[0*4+2]*inMatrix[2*4+3]*inMatrix[3*4+1] - inMatrix[0*4+3]*inMatrix[2*4+1]*inMatrix[3*4+2]);
	resultMatrix[0*4+2] = invDeterminant * (inMatrix[0*4+1]*inMatrix[1*4+2]*inMatrix[3*4+3] + inMatrix[0*4+2]*inMatrix[1*4+3]*inMatrix[3*4+1] + inMatrix[0*4+3]*inMatrix[1*4+1]*inMatrix[3*4+2] - inMatrix[0*4+1]*inMatrix[1*4+3]*inMatrix[3*4+2] - inMatrix[0*4+2]*inMatrix[1*4+1]*inMatrix[3*4+3] - inMatrix[0*4+3]*inMatrix[1*4+2]*inMatrix[3*4+1]);
	resultMatrix[0*4+3] = invDeterminant * (inMatrix[0*4+1]*inMatrix[1*4+3]*inMatrix[2*4+2] + inMatrix[0*4+2]*inMatrix[1*4+1]*inMatrix[2*4+3] + inMatrix[0*4+3]*inMatrix[1*4+2]*inMatrix[2*4+1] - inMatrix[0*4+1]*inMatrix[1*4+2]*inMatrix[2*4+3] - inMatrix[0*4+2]*inMatrix[1*4+3]*inMatrix[2*4+1] - inMatrix[0*4+3]*inMatrix[1*4+1]*inMatrix[2*4+2]);
	resultMatrix[1*4+0] = invDeterminant * (inMatrix[1*4+0]*inMatrix[2*4+3]*inMatrix[3*4+2] + inMatrix[1*4+2]*inMatrix[2*4+0]*inMatrix[3*4+3] + inMatrix[1*4+3]*inMatrix[2*4+2]*inMatrix[3*4+0] - inMatrix[1*4+0]*inMatrix[2*4+2]*inMatrix[3*4+3] - inMatrix[1*4+2]*inMatrix[2*4+3]*inMatrix[3*4+0] - inMatrix[1*4+3]*inMatrix[2*4+0]*inMatrix[3*4+2]);
	resultMatrix[1*4+1] = invDeterminant * (inMatrix[0*4+0]*inMatrix[2*4+2]*inMatrix[3*4+3] + inMatrix[0*4+2]*inMatrix[2*4+3]*inMatrix[3*4+0] + inMatrix[0*4+3]*inMatrix[2*4+0]*inMatrix[3*4+2] - inMatrix[0*4+0]*inMatrix[2*4+3]*inMatrix[3*4+2] - inMatrix[0*4+2]*inMatrix[2*4+0]*inMatrix[3*4+3] - inMatrix[0*4+3]*inMatrix[2*4+2]*inMatrix[3*4+0]);
	resultMatrix[1*4+2] = invDeterminant * (inMatrix[0*4+0]*inMatrix[1*4+3]*inMatrix[3*4+2] + inMatrix[0*4+2]*inMatrix[1*4+0]*inMatrix[3*4+3] + inMatrix[0*4+3]*inMatrix[1*4+2]*inMatrix[3*4+0] - inMatrix[0*4+0]*inMatrix[1*4+2]*inMatrix[3*4+3] - inMatrix[0*4+2]*inMatrix[1*4+3]*inMatrix[3*4+0] - inMatrix[0*4+3]*inMatrix[1*4+0]*inMatrix[3*4+2]);
	resultMatrix[1*4+3] = invDeterminant * (inMatrix[0*4+0]*inMatrix[1*4+2]*inMatrix[2*4+3] + inMatrix[0*4+2]*inMatrix[1*4+3]*inMatrix[2*4+0] + inMatrix[0*4+3]*inMatrix[1*4+0]*inMatrix[2*4+2] - inMatrix[0*4+0]*inMatrix[1*4+3]*inMatrix[2*4+2] - inMatrix[0*4+2]*inMatrix[1*4+0]*inMatrix[2*4+3] - inMatrix[0*4+3]*inMatrix[1*4+2]*inMatrix[2*4+0]);
	resultMatrix[2*4+0] = invDeterminant * (inMatrix[1*4+0]*inMatrix[2*4+1]*inMatrix[3*4+3] + inMatrix[1*4+1]*inMatrix[2*4+3]*inMatrix[3*4+0] + inMatrix[1*4+3]*inMatrix[2*4+0]*inMatrix[3*4+1] - inMatrix[1*4+0]*inMatrix[2*4+3]*inMatrix[3*4+1] - inMatrix[1*4+1]*inMatrix[2*4+0]*inMatrix[3*4+3] - inMatrix[1*4+3]*inMatrix[2*4+1]*inMatrix[3*4+0]);
	resultMatrix[2*4+1] = invDeterminant * (inMatrix[0*4+0]*inMatrix[2*4+3]*inMatrix[3*4+1] + inMatrix[0*4+1]*inMatrix[2*4+0]*inMatrix[3*4+3] + inMatrix[0*4+3]*inMatrix[2*4+1]*inMatrix[3*4+0] - inMatrix[0*4+0]*inMatrix[2*4+1]*inMatrix[3*4+3] - inMatrix[0*4+1]*inMatrix[2*4+3]*inMatrix[3*4+0] - inMatrix[0*4+3]*inMatrix[2*4+0]*inMatrix[3*4+1]);
	resultMatrix[2*4+2] = invDeterminant * (inMatrix[0*4+0]*inMatrix[1*4+1]*inMatrix[3*4+3] + inMatrix[0*4+1]*inMatrix[1*4+3]*inMatrix[3*4+0] + inMatrix[0*4+3]*inMatrix[1*4+0]*inMatrix[3*4+1] - inMatrix[0*4+0]*inMatrix[1*4+3]*inMatrix[3*4+1] - inMatrix[0*4+1]*inMatrix[1*4+0]*inMatrix[3*4+3] - inMatrix[0*4+3]*inMatrix[1*4+1]*inMatrix[3*4+0]);
	resultMatrix[2*4+3] = invDeterminant * (inMatrix[0*4+0]*inMatrix[1*4+3]*inMatrix[2*4+1] + inMatrix[0*4+1]*inMatrix[1*4+0]*inMatrix[2*4+3] + inMatrix[0*4+3]*inMatrix[1*4+1]*inMatrix[2*4+0] - inMatrix[0*4+0]*inMatrix[1*4+1]*inMatrix[2*4+3] - inMatrix[0*4+1]*inMatrix[1*4+3]*inMatrix[2*4+0] - inMatrix[0*4+3]*inMatrix[1*4+0]*inMatrix[2*4+1]);
	resultMatrix[3*4+0] = invDeterminant * (inMatrix[1*4+0]*inMatrix[2*4+2]*inMatrix[3*4+1] + inMatrix[1*4+1]*inMatrix[2*4+0]*inMatrix[3*4+2] + inMatrix[1*4+2]*inMatrix[2*4+1]*inMatrix[3*4+0] - inMatrix[1*4+0]*inMatrix[2*4+1]*inMatrix[3*4+2] - inMatrix[1*4+1]*inMatrix[2*4+2]*inMatrix[3*4+0] - inMatrix[1*4+2]*inMatrix[2*4+0]*inMatrix[3*4+1]);
	resultMatrix[3*4+1] = invDeterminant * (inMatrix[0*4+0]*inMatrix[2*4+1]*inMatrix[3*4+2] + inMatrix[0*4+1]*inMatrix[2*4+2]*inMatrix[3*4+0] + inMatrix[0*4+2]*inMatrix[2*4+0]*inMatrix[3*4+1] - inMatrix[0*4+0]*inMatrix[2*4+2]*inMatrix[3*4+1] - inMatrix[0*4+1]*inMatrix[2*4+0]*inMatrix[3*4+2] - inMatrix[0*4+2]*inMatrix[2*4+1]*inMatrix[3*4+0]);
	resultMatrix[3*4+2] = invDeterminant * (inMatrix[0*4+0]*inMatrix[1*4+2]*inMatrix[3*4+1] + inMatrix[0*4+1]*inMatrix[1*4+0]*inMatrix[3*4+2] + inMatrix[0*4+2]*inMatrix[1*4+1]*inMatrix[3*4+0] - inMatrix[0*4+0]*inMatrix[1*4+1]*inMatrix[3*4+2] - inMatrix[0*4+1]*inMatrix[1*4+2]*inMatrix[3*4+0] - inMatrix[0*4+2]*inMatrix[1*4+0]*inMatrix[3*4+1]);
	resultMatrix[3*4+3] = invDeterminant * (inMatrix[0*4+0]*inMatrix[1*4+1]*inMatrix[2*4+2] + inMatrix[0*4+1]*inMatrix[1*4+2]*inMatrix[2*4+0] + inMatrix[0*4+2]*inMatrix[1*4+0]*inMatrix[2*4+1] - inMatrix[0*4+0]*inMatrix[1*4+2]*inMatrix[2*4+1] - inMatrix[0*4+1]*inMatrix[1*4+0]*inMatrix[2*4+2] - inMatrix[0*4+2]*inMatrix[1*4+1]*inMatrix[2*4+0]);

}