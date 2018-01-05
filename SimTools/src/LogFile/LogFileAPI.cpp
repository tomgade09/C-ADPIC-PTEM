#include "LogFile\LogFileAPI.h"

DLLEXPORT void writeLogFileEntryAPI(LogFile* log, const char* logMessage) {
	log->writeLogFileEntry(logMessage); }

DLLEXPORT void writeTimeDiffFromNowAPI(LogFile* log, int startTSind, const char* nowLabel) {
	log->writeTimeDiffFromNow(startTSind, nowLabel); }

DLLEXPORT void writeTimeDiffAPI(LogFile* log, int startTSind, int endTSind) {
	log->writeTimeDiff(startTSind, endTSind); }