#include "LogFile\LogFile.h"
#include <iomanip>

void LogFile::writeLogFileEntry(std::string logMessage)
{//[Time (ms from start) - 19 chars tot, 15 chars for numbers] | Log Data - 20 chars | Log Message - unlimited chars
	std::string writeTxt;
	std::stringstream ss;

	ss << std::setprecision(10) << std::setw(11) << static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timeStructs_m[0]->tp).count()) / 1000;
	writeTxt = "[ " + ss.str() + " ] : " + logMessage + "\n"; //time
	
	fileIO::writeTxtFile(logFileName_m.c_str(), writeTxt.c_str());
}

void LogFile::writeErrorEntry(std::string functionName, std::string logMessage, std::vector<std::string> args)
{
	std::string error{ "ERROR: " + functionName + "(" };
	for (size_t arg = 0; arg < args.size(); arg++)
		error += args.at(arg) + ((arg == args.size()-1) ? "" : ", ");
	error += ")";

	writeLogFileEntry(error);

	logMessage = "                  -> " + logMessage + "\n";

	fileIO::writeTxtFile(logFileName_m.c_str(), logMessage.c_str());
}

void LogFile::createTimeStruct(std::string label)
{
	timeStruct* tS = new timeStruct;
	tS->label = label;
	tS->tp = std::chrono::steady_clock::now();
	timeStructs_m.push_back(tS);
}

//writeTimeDiff plus overloads
void LogFile::writeTimeDiff(timeStruct* startTS, timeStruct* endTS)
{
	std::string logMessage;
	logMessage = "LogFile::writeTimeDiff:  " + startTS->label + "  TO  " + endTS->label + ": " \
		+ std::to_string(static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(endTS->tp - startTS->tp).count()) / 1000);

	writeLogFileEntry(logMessage);
}

void LogFile::writeTimeDiff(int startTSind, timeStruct* endTS)
{
	writeTimeDiff(timeStructs_m[startTSind], endTS);
}

void LogFile::writeTimeDiff(int startTSind, int endTSind)
{
	writeTimeDiff(timeStructs_m[startTSind], timeStructs_m[endTSind] );
}


//writeTimeDiffFromNow plus overloads
void LogFile::writeTimeDiffFromNow(timeStruct* startTS, std::string nowLabel)
{
	timeStruct* ts = new timeStruct;
	ts->label = nowLabel;
	ts->tp = std::chrono::steady_clock::now();

	writeTimeDiff(startTS, ts);
}

void LogFile::writeTimeDiffFromNow(int startTSind, std::string nowLabel)
{
	writeTimeDiffFromNow(timeStructs_m[startTSind], nowLabel);
}


void LogFile::printTimeNowFromLastTS()
{
	std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timeStructs_m[timeStructs_m.size()-1]->tp).count()) / 1000000;
}