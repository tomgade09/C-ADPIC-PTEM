#ifndef SA_BINARYFILETOOLS_H
#define SA_BINARYFILETOOLS_H

#include "FileIO\fileIO.h"

std::vector<std::string> constCharToStrVec(const char* str, const char delim=',');
void stringPadder(std::string& in, int totalStrLen, int indEraseFrom = 1);
double*** form3Darray(int d1, int d2, int d3);
void delete3Darray(double*** array3D, int d1);
std::vector<std::vector<std::vector<double>>> form3DvectorArray(int numPartTypes, int numAttr, int numParts);
std::vector<std::vector<double>> form2DvectorArray(int numAttr, int numParts);

#endif