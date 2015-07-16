#include "Cylinder.h"
#include "globals.h"


//Initialize the date for a cylinder
//baseRadius - radius of bottom section
//topRadius - radius of top section
//height - height of cylinder
//horizSlices - horizontal slices (portions of a circle)
//vertStacks - vertical slices (determines number of circles used)
//color - the color used
void Cylinder::init(double baseRadius, double topRadius, double height, int horizSlices, int vertStacks, vec3 color1, vec3 color2)
{
		
	vertices.clear();
	indices.clear();

	int indexCount = 0; //Current vertex index counter

	
	double increment1 = 2 * 3.14159 / horizSlices;
	double increment2 = height / vertStacks;

	double radiusIncrement = (topRadius - baseRadius) / vertStacks;
	

	//Iterate through all angles within a circle
	for (double angle = 0; angle <= 2 * 3.14159 + 0.0001; angle += increment1)
	{
		//Add one slice on the bottom of the cylinder
		double x1 = cos(angle);
		double z1 = sin(angle);
		double x2 = cos(angle + increment1);
		double z2 = sin(angle + increment1);

		vec3 normalA(0, -height / 2, 0);
		normalA = normalize(normalA);
		VertexB vertexA = {0, -height / 2, 0, 1, normalA[0], normalA[1], normalA[2], 0, color1[0], color1[1], color1[2], 1.0f};
		vec3 normalB(x1*baseRadius, -height / 2, z1*baseRadius);
		normalB = normalize(normalB);
		VertexB vertexB = {x1*baseRadius, -height / 2, z1*baseRadius, 1, normalB[0], normalB[1], normalB[2], 0, color2[0], color2[1], color2[2], 1.0f};
		vec3 normalC(x2*baseRadius, -height / 2, z2*baseRadius);
		normalC = normalize(normalC);
		VertexB vertexC = {x2*baseRadius, -height / 2, z2*baseRadius, 1, normalC[0], normalC[1], normalC[2], 0, color2[0], color2[1], color2[2], 1.0f};

		vertices.push_back(vertexA);
		vertices.push_back(vertexB);
		vertices.push_back(vertexC);

		indices.push_back(indexCount);
		indices.push_back(indexCount+1);
		indices.push_back(indexCount+2);
		indexCount += 3;


		//For the current circle sice, add it to the mesh at each height increment
		for (double h = -height/2, r = baseRadius; h <= height / 2 - increment2 + 0.0001; h += increment2, r += radiusIncrement)
		{
			//Struct initialization syntax to make code more compact
			vec3 normal1(x1 * r, h, z1 * r);
			normal1 = normalize(normal1);
			VertexB vertex1 = {x1 * r, h, z1 * r, 1, normal1[0], normal1[1], normal1[2], 0, color1[0], color1[1], color1[2], 1.0f};
			vec3 normal2(x1 * (r + radiusIncrement), h + increment2, z1 * (r + radiusIncrement));
			normal2 = normalize(normal2);
			VertexB vertex2 = {x1 * (r + radiusIncrement), h + increment2, z1 * (r + radiusIncrement), 1, normal2[0], normal2[1], normal2[2], 0, color2[0], color2[1], color2[2], 1.0f};
			vec3 normal3(x2 * (r + radiusIncrement), h + increment2, z2 * (r + radiusIncrement));
			normal3 = normalize(normal3);
			VertexB vertex3 = {x2 * (r + radiusIncrement), h + increment2, z2 * (r + radiusIncrement), 1, normal3[0], normal3[1], normal3[2], 0, color2[0], color2[1], color2[2], 1.0f};
			vec3 normal4(x2 * r, h, z2 * r);
			normal4 = normalize(normal4);
			VertexB vertex4 = {x2 * r, h, z2 * r, 1, normal4[0], normal4[1], normal4[2], 0, color1[0], color1[1], color1[2], 1.0f};


			vertices.push_back(vertex1);
			vertices.push_back(vertex2);
			vertices.push_back(vertex3);
			vertices.push_back(vertex4);

			indices.push_back(indexCount);
			indices.push_back(indexCount+1);
			indices.push_back(indexCount+2);
			indices.push_back(indexCount+0);
			indices.push_back(indexCount+3);
			indices.push_back(indexCount+2);

			indexCount += 4;

		}

		//Add the current circle slice to the top of the cylinder
		vec3 normalD(0, height / 2, 0);
		normalD = normalize(normalD);
		VertexB vertexD = {0, height / 2, 0, 1, normalD[0], normalD[1], normalD[2], 0, color1[0], color1[1], color1[2], 1.0f};
		vec3 normalE(x1*topRadius, height / 2, z1*topRadius);
		normalE = normalize(normalE);
		VertexB vertexE = {x1*topRadius, height / 2, z1*topRadius, 1, normalE[0], normalE[1], normalE[2], 0, color2[0], color2[1], color2[2], 1.0f};
		vec3 normalF(x2*topRadius, height / 2, z2*topRadius);
		normalF = normalize(normalF);
		VertexB vertexF = {x2*topRadius, height / 2, z2*topRadius, 1, normalF[0], normalF[1], normalF[2], 0, color2[0], color2[1], color2[2], 1.0f};

		vertices.push_back(vertexD);
		vertices.push_back(vertexE);
		vertices.push_back(vertexF);

		indices.push_back(indexCount);
		indices.push_back(indexCount+1);
		indices.push_back(indexCount+2);
		indexCount += 3;
	}
	
}

