#ifndef LOGFILEAPI_H
#define LOGFILEAPI_H

#include "dlldefines.h"
#include "LogFile/LogFile.h"

DLLEXP_EXTC void writeLogFileEntryAPI(LogFile* log, const char* logMessage);
DLLEXP_EXTC void writeTimeDiffFromNowAPI(LogFile* log, int startTSind, const char* nowLabel);
DLLEXP_EXTC void writeTimeDiffAPI(LogFile* log, int startTSind, int endTSind);

#endif