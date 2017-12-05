#ifndef SA_BINARYFILETOOLS_H
#define SA_BINARYFILETOOLS_H

#include "FileIO\fileIO.h"

void saveParticleAttributeToDisk(double* arrayToSave, int length, const char* foldername, const char* name);
void loadFileIntoParticleAttribute(double* arrayToLoadInto, int length, const char* foldername, const char* name);

#endif