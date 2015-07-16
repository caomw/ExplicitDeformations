#include "Cube.h"
#include "globals.h"


//Initialize the data for a cube
//width/height/depth - size parameters of cube
//color - the color applied
void Cube::init(double width, double height, double depth, vec3 color1, vec3 color2)
{

	
	vertices.clear();
	indices.clear();

	//Half of each dimension so that it is centered on (0, 0, 0)
	double halfWidth = width / 2;
	double halfHeight = height / 2;
	double halfDepth = depth / 2;

	vec3 normal0(-halfWidth, -halfHeight, halfDepth);
	normal0 = normalize(normal0);
	VertexB vertex0 = {-halfWidth, -halfHeight, halfDepth, 1, normal0[0], normal0[1], normal0[2], 0, color1[0], color1[1], color1[2], 1.0f}; //Front lower left
	vec3 normal1(-halfWidth, halfHeight, halfDepth);
	normal1 = normalize(normal1);
	VertexB vertex1 = {-halfWidth, halfHeight, halfDepth, 1, normal1[0], normal1[1], normal1[2], 0, color2[0], color2[1], color2[2], 1.0f}; //Front upper left
	vec3 normal2(halfWidth, halfHeight, halfDepth);
	normal2 = normalize(normal2);
	VertexB vertex2 = {halfWidth, halfHeight, halfDepth, 1, normal2[0], normal2[1], normal2[2], 0, color1[0], color1[1], color1[2], 1.0f}; //Front top right
	vec3 normal3(halfWidth, -halfHeight, halfDepth);
	normal3 = normalize(normal3);
	VertexB vertex3 = {halfWidth, -halfHeight, halfDepth, 1, normal3[0], normal3[1], normal3[2], 0, color2[0], color2[1], color2[2], 1.0f}; //Front lower right


	vec3 normal4(-halfWidth, -halfHeight, -halfDepth);
	normal4 = normalize(normal4);

	VertexB vertex4 = {-halfWidth, -halfHeight, -halfDepth, 1, normal4[0], normal4[1], normal4[2], 0, color1[0], color1[1], color1[2], 1.0f}; //Back lower left
	
	vec3 normal5(-halfWidth, halfHeight, -halfDepth);
	normal5 = normalize(normal5);
	
	VertexB vertex5 = {-halfWidth, halfHeight, -halfDepth, 1, normal5[0], normal5[1], normal5[2], 0, color2[0], color2[1], color2[2], 1.0f}; //Back upper left
	
	vec3 normal6(halfWidth, halfHeight, -halfDepth);
	normal6 = normalize(normal6);
	
	VertexB vertex6 = {halfWidth, halfHeight, -halfDepth, 1,  normal6[0], normal6[1], normal6[2], 0, color1[0], color1[1], color1[2], 1.0f}; //Back top right
	
	vec3 normal7(halfWidth, -halfHeight, -halfDepth);
	normal7 = normalize(normal7);
	
	VertexB vertex7 = {halfWidth, -halfHeight, -halfDepth, 1,  normal7[0], normal7[1], normal7[2], 0, color2[0], color2[1], color2[2], 1.0f}; //Back lower right
	
	vertices.push_back(vertex0);
	vertices.push_back(vertex1);
	vertices.push_back(vertex2);
	vertices.push_back(vertex3);
	vertices.push_back(vertex4);
	vertices.push_back(vertex5);
	vertices.push_back(vertex6);
	vertices.push_back(vertex7);

	//Front quad
	indices.push_back(0);
	indices.push_back(1);
	indices.push_back(2);
	indices.push_back(2);
	indices.push_back(3);
	indices.push_back(0);
	
	//Left Quad
	indices.push_back(4);
	indices.push_back(5);
	indices.push_back(1);
	indices.push_back(1);
	indices.push_back(0);
	indices.push_back(4);

	//Back quad
	indices.push_back(4);
	indices.push_back(5);
	indices.push_back(6);
	indices.push_back(6);
	indices.push_back(7);
	indices.push_back(4);

	//Right quad
	indices.push_back(3);
	indices.push_back(2);
	indices.push_back(6);
	indices.push_back(6);
	indices.push_back(7);
	indices.push_back(3);
	

	//Bottom Quad
	indices.push_back(4);
	indices.push_back(0);
	indices.push_back(3);
	indices.push_back(3);
	indices.push_back(7);
	indices.push_back(4);
	
	//Top Quad
	indices.push_back(1);
	indices.push_back(5);
	indices.push_back(6);
	indices.push_back(6);
	indices.push_back(2);
	indices.push_back(1);

}


