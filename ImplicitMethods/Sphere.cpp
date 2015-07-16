#include "Sphere.h"
#include "globals.h"

//Initialize sphere mesh
//radius - radius of sphere
//HorizSlices - Number of horizontal slices to generate
//VertStacks - Number of vertical stacks to generate
//color - the color to apply
void Sphere::init(double radius, int horizSlices, int vertStacks, vec3 color1, vec3 color2)
{	
	vertices.clear();
	indices.clear();

	int indexCount = 0;
	
	double increment1 = 2 * 3.14159 / horizSlices;
	double increment2 = 3.14159 / vertStacks;

	for (double theta = 0; theta <= 2 * 3.14159 + 0.0001; theta += increment1)
	{
		for (double phi = 0; phi <= 3.14159 + 0.0001; phi += increment2)
		{
			vec4 transformedPoint1 = getSphericalPoint(radius, theta, phi);
			vec3 normal1(transformedPoint1[0], transformedPoint1[1], transformedPoint1[2]);
			normal1 = normalize(normal1);
			vec4 transformedPoint2 = getSphericalPoint(radius, theta + increment1, phi);	
			vec3 normal2(transformedPoint2[0], transformedPoint2[1], transformedPoint2[2]);
			normal2 = normalize(normal2);
			vec4 transformedPoint3 = getSphericalPoint(radius, theta, phi + increment2);
			vec3 normal3(transformedPoint3[0], transformedPoint3[1], transformedPoint3[2]);
			normal3 = normalize(normal3);
			vec4 transformedPoint4 = getSphericalPoint(radius, theta + increment1, phi + increment2);
			vec3 normal4(transformedPoint4[0], transformedPoint4[1], transformedPoint4[2]);
			normal4 = normalize(normal4);

			//Struct initialization syntax to make code more compact
			VertexB vertex1 = {transformedPoint1[0], transformedPoint1[1], transformedPoint1[2], transformedPoint1[3], normal1[0], normal1[1], normal1[2], 0, color1[0], color1[1], color1[2], 1.0f};
			VertexB vertex2 = {transformedPoint2[0], transformedPoint2[1], transformedPoint2[2], transformedPoint2[3], normal2[0], normal2[1], normal2[2], 0, color2[0], color2[1], color2[2], 1.0f};
			VertexB vertex3 = {transformedPoint3[0], transformedPoint3[1], transformedPoint3[2], transformedPoint3[3], normal3[0], normal3[1], normal3[2], 0, color2[0], color2[1], color2[2], 1.0f};
			VertexB vertex4 = {transformedPoint4[0], transformedPoint4[1], transformedPoint4[2], transformedPoint4[3], normal4[0], normal4[1], normal4[2], 0, color1[0], color1[1], color1[2], 1.0f};


			vertices.push_back(vertex1);
			vertices.push_back(vertex2);
			vertices.push_back(vertex3);
			vertices.push_back(vertex4);

			indices.push_back(indexCount);
			indices.push_back(indexCount+1);
			indices.push_back(indexCount+2);
			indices.push_back(indexCount+3);
			indices.push_back(indexCount+1);
			indices.push_back(indexCount+2);

			indexCount += 4;

		}
	}
	
}


//This function returns the cartesian coordinates for a given spherical coordinate
//Radius - radius of sphere
//theta, phi - spherical coordinates
//Based on equations at http://mathworld.wolfram.com/SphericalCoordinates.html
vec4 Sphere::getSphericalPoint(double radius, double theta, double phi)
{
	vec4 sphericalPoint(radius * cos(theta) * sin(phi), radius * sin(theta) * sin(phi), radius * cos(phi), 1);
	vec4 transformedPoint = sphericalPoint;
	return transformedPoint;

}
