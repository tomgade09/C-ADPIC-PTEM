#ifndef LOGFILE_H
#define LOGFILE_H

#include <chrono>
#include <string>
#include <vector>
#include <iomanip>
#include "FileIO\fileIO.h"

struct timeStruct
{
	std::string label;
	std::chrono::steady_clock::time_point tp;
};

class LogFile
{
private:
	std::string logFileName_m;
	bool overwrite_m;
	std::vector<timeStruct*> timeStructs_m;

public:
	LogFile(std::string logFileName, int timeStructToReserve, bool overwrite=false) : logFileName_m{ logFileName }, overwrite_m{ overwrite }
	{ 
		timeStructs_m.reserve(timeStructToReserve);
		std::string logHeader{ "   [ Time  (ms) ]   |       Log Data       | Log Message\n\n" }; //do I want to add the time, other attributes to file???
		logHeader += "[ 0 ]               | LogFile class        | Log file class created, file created on disk, first entry written, first time point recorded.\n";
		fileIO::writeTxtFile(logFileName_m.c_str(), logHeader.c_str(), overwrite_m);
		createTimeStruct("LogFile class constructor"); //index 0 of timeStructs_m
	}

	~LogFile()
	{
		for (int iii = 0; iii < timeStructs_m.size(); iii++)
			delete timeStructs_m[iii];
	}

	void createTimeStruct(std::string label);
	void writeLogFileEntry(std::string logData, std::string logMessage);
	
	void writeTimeDiff(timeStruct* startTS, timeStruct* endTS);
	void writeTimeDiff(int startTSind, timeStruct* endTS);
	void writeTimeDiff(int startTSind, int endTSind);
	
	void writeTimeDiffFromNow(timeStruct* startTS, std::string nowLabel);
	void writeTimeDiffFromNow(int startTSind, std::string nowLabel);
	
	void printTimeNowFromLastTS();
};

#endif