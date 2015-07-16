#ifndef GLOBALS_H
#define GLOBALS_H

//Global variables
extern double width; //Width of one cell
extern double height; //Depth of one cell
extern GLuint programObject; //GLSL shading program object
extern bool fillMode; //True if filling polygons; false if using lines
extern bool ambientMode; //True if using high ambient mode (for debugging); false if normal ambience

#endif