#pragma once

#include "ParticleSystem.h"

//Particle System class - Stanford method paper
//It builds a grid of particles (represented by Particle objects)
//These are then held together by springs using hooke's laws.  One edge corresponds to each spring.
//Hooke's law is used to represent the springs, with both a spring component and a damping comopnent
//Implicit methods are used for the integration.  This requires a solving system but allows much larger spring constants without instability
//and potentially allows more efficient implementation
//A gravity component is also present
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
