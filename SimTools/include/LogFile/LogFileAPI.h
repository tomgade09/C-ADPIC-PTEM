#ifndef LOGFILEAPI_H
#define LOGFILEAPI_H

#include "LogFile\LogFile.h"

DLLEXPORT void writeLogFileEntryAPI(LogFile* log, const char* logData, const char* logMessage);
DLLEXPORT void writeTimeDiffFromNowAPI(LogFile* log, int startTSind, const char* nowLabel);
DLLEXPORT void writeTimeDiffAPI(LogFile* log, int startTSind, int endTSind);

#endif