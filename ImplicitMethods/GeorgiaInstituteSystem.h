#pragma once

#include "ParticleSystem.h"

//Particle System class - Deformation Method #2
//
//Based on the paper at http://graphics.berkeley.edu/papers/Obrien-GMA-1999-08/Obrien-GMA-1999-08.pdf – Graphical Modeling and Animation of Brittle Fracture
//By James O'Brien and Jessica Hodgkins
class GeorgiaInstituteSystem : public ParticleSystem
{
	public:
		GeorgiaInstituteSystem(Vertex * vertexList, int vertexCount, int * tetraList, int tetraCount, Logger * logger);
		~GeorgiaInstituteSystem();
		void doUpdate(double elapsedSeconds);
	private:
		double * m;
		double * beta;

};
