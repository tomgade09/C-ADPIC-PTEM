#ifndef LOGFILEAPI_H
#define LOGFILEAPI_H

#include "dllexport.h"
#include "LogFile\LogFile.h"

DLLEXPORT void writeLogFileEntryAPI(LogFile* log, const char* logMessage);
DLLEXPORT void writeTimeDiffFromNowAPI(LogFile* log, int startTSind, const char* nowLabel);
DLLEXPORT void writeTimeDiffAPI(LogFile* log, int startTSind, int endTSind);

#endif