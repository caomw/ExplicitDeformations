#pragma once

#include "ParticleSystem.h"

//Particle System class - Deformation Method #1
//
//Based on the paper at http://www.math.ucla.edu/~jteran/papers/TSNF03.pdf – Finite Volume Methods for the Simulation of Skeletal Muscle
//By R. Fedkiw et. al
class StanfordSystem : public ParticleSystem
{
	public:
	StanfordSystem(Vertex * vertexList, int vertexCount, int * tetraList, int tetraCount, Logger * logger);
	~StanfordSystem();
	double * crossProductSums;
	double * invDm;
	void doUpdate(double elapsedSeconds);
};
