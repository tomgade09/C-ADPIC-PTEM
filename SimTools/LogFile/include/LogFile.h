#ifndef LOGFILE_H
#define LOGFILE_H

#include <chrono>
#include <string>
#include <vector>
#include "include\fileIO.h"

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
	std::vector<timeStruct> timeStructs_m;

public:
	LogFile(std::string logFileName, int timeStructToReserve, bool overwrite=false) : logFileName_m{ logFileName }, overwrite_m{ overwrite }
	{ 
		timeStructs_m.reserve(timeStructToReserve);
		std::string logHeader{ " [ Time (ms) ]  |    Log  Data    | Log Message\n\n" }; //do I want to add the time, other attributes to file???
		logHeader += "[ 0 ]           | Create Log File | Log file class created, file created on disk, first entry written, first time point recorded.\n";
		fileIO::writeTxtFile(logFileName_m.c_str(), logHeader.c_str(), overwrite_m);
		createTimeStruct("LogFile class constructor"); //index 0 of timeStructs_m
	}

	~LogFile()
	{

	}

	void createTimeStruct(std::string label);
	void writeLogFileEntry(std::string logData, std::string logMessage);
	
	void writeLogTimeDiff(timeStruct startTS, timeStruct endTS);
	void writeLogTimeDiff(int startTSind, timeStruct endTS);
	void writeLogTimeDiff(int startTSind, int endTSind);
	
	void writeLogTimeDiffFromNow(timeStruct startTS, std::string nowLabel);
	void writeLogTimeDiffFromNow(int startTSind, std::string nowLabel);
	
	
	
	
	
	
	void printTimeNowFromTimeStruct(timeStruct* tS, std::string label);
	void printTimeNowFromTSJustMS(timeStruct* tS);
	void printTimeDiffBtwTwoTimeStructs(timeStruct* startTS, timeStruct* endTS);
	void printTimeDiffJustTimeMS(timeStruct* startTS, timeStruct* endTS);
};

#endif