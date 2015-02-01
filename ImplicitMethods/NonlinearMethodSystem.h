#pragma once

#include "ParticleSystem.h"

//Particle System class - Deformation Method #3
//
//Based on the paper at http://www-ljk.imag.fr/Publications/Basilic/com.lmc.publi.PUBLI_Article@11f6a0378d9_18c74/tensile.pdf – Simple, yet Accurate Nonlinear Tensile Stiffness
//Pascal Volino et. al
class NonlinearMethodSystem : public ParticleSystem
{
	public:
	NonlinearMethodSystem(Vertex * vertexList, int vertexCount, int * tetraList, int tetraCount, Logger * logger);
	~NonlinearMethodSystem();
	void doUpdate(double elapsedSeconds);

	double * ruWeights;
	double * rvWeights;
	double * rwWeights;
};
