#ifndef SA_BINARYFILETOOLS_H
#define SA_BINARYFILETOOLS_H

#include "FileIO\fileIO.h"

void saveParticleAttributeToDisk(double* arrayToSave, int length, const char* foldername, const char* name);
void loadFileIntoParticleAttribute(double* arrayToLoadInto, int length, const char* foldername, const char* name);
void stringPadder(std::string& in, int totalStrLen, int indEraseFrom = 1);
double*** form3Darray(int d1, int d2, int d3);
void delete3Darray(double*** array3D, int d1);

#endif