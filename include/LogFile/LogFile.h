#ifndef LOGFILE_H
#define LOGFILE_H

#include <string>
#include <vector>
#include <memory>
#include <chrono>

struct timeStruct
{
	std::string label;
	std::chrono::steady_clock::time_point tp;

	timeStruct(std::string lbl, std::chrono::steady_clock::time_point time) :
		label{ lbl }, tp{ time } {}
};

class LogFile
{
private:
	std::string logFileName_m;
	bool overwrite_m;
	std::vector<std::unique_ptr<timeStruct>> timeStructs_m;

public:
	LogFile(std::string logFileName, int timeStructToReserve, bool overwrite = false);

	~LogFile() {}

	void createTimeStruct(std::string label);
	void writeLogFileEntry(std::string logMessage);
	
	void writeTimeDiff(timeStruct* startTS, timeStruct* endTS);
	void writeTimeDiff(int startTSind, timeStruct* endTS);
	void writeTimeDiff(int startTSind, int endTSind);
	
	void writeTimeDiffFromNow(timeStruct* startTS, std::string nowLabel);
	void writeTimeDiffFromNow(int startTSind, std::string nowLabel);
	
	void printTimeNowFromLastTS();
	void printTimeNowFromFirstTS();
};

#endif