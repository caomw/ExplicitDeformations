#pragma once

const int DIMENSION = 3;	//Number of dimensions - 3 means 3D

//This struct represents one vertex in a tetrahedral mesh
struct Vertex
{
	Vertex()
	{
		position[DIMENSION] = 1;
		vertexNormal[DIMENSION] = 1;
		color[0] = 0.0;
		color[1] = 1.0;
		color[2] = 0.0;
		color[3] = 1.0;
		velocity[DIMENSION] = 1;
	}
	Vertex (float p1, float p2, float p3, float p4, float p5, float p6, float p7, float p8, float p9, float p10, float p11, float p12, float p13, float p14, float p15, float p16, int p17)
	{
		position[0] = p1;
		position[1] = p2;
		position[2] = p3;
		position[3] = p4;
		vertexNormal[0] = p5;
		vertexNormal[1] = p6;
		vertexNormal[2] = p7;
		vertexNormal[3] = p8;
		color[0] = p9;
		color[1] = p10;
		color[2] = p11;
		color[3] = p12;
		velocity[0] = p13;
		velocity[1] = p14;
		velocity[2] = p15;
		velocity[3] = p16;
		triangleCount = p17;
	}

	float position [DIMENSION + 1];	//Position vector in 3d space
	float vertexNormal [DIMENSION + 1];		//Normal vector in 3d space (this one is a VERTEX AVERAGE used in lighting)
	float color[4];
	float velocity [DIMENSION + 1];	//Velocity vector in 3d space	
	int triangleCount;				//Number of triangles in which the vertex for this particle is contained (helps calculate normals)
};

