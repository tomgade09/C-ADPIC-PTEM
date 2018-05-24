#include "API/LogFileAPI.h"

DLLEXP_EXTC void writeLogFileEntryAPI(LogFile* log, const char* logMessage) {
	log->writeLogFileEntry(logMessage); }

DLLEXP_EXTC void writeTimeDiffFromNowAPI(LogFile* log, int startTSind, const char* nowLabel) {
	log->writeTimeDiffFromNow(startTSind, nowLabel); }

DLLEXP_EXTC void writeTimeDiffAPI(LogFile* log, int startTSind, int endTSind) {
	log->writeTimeDiff(startTSind, endTSind); }