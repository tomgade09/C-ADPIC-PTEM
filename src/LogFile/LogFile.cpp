#include "LogFile\LogFile.h"
#include <iomanip>

void LogFile::writeLogFileEntry(std::string logMessage)
{//[Time (ms from start) - 19 chars tot, 15 chars for numbers] | Log Data - 20 chars | Log Message - unlimited chars
	std::string writeTxt;
	std::stringstream ss;

	ss << std::setprecision(10) << std::setw(11) << static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timeStructs_m.at(0)->tp).count()) / 1000;
	std::string ssstr{ ss.str() };
	while (ssstr.length() > 11)
		ssstr.erase(std::ios::end);
	writeTxt = "[ " + ssstr + " ] : " + logMessage + "\n"; //time

	FILE_RDWR_EXCEP_CHECK(fileIO::writeTxtFile(writeTxt, logFileName_m));
}

void LogFile::createTimeStruct(std::string label)
{
	timeStructs_m.push_back(std::make_unique<timeStruct>(label, std::chrono::steady_clock::now()));
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
	writeTimeDiff(timeStructs_m.at(startTSind).get(), endTS);
}

void LogFile::writeTimeDiff(int startTSind, int endTSind)
{
	writeTimeDiff(timeStructs_m.at(startTSind).get(), timeStructs_m.at(endTSind).get() );
}


//writeTimeDiffFromNow plus overloads
void LogFile::writeTimeDiffFromNow(timeStruct* startTS, std::string nowLabel)
{
	std::unique_ptr<timeStruct> ts{ std::make_unique<timeStruct>(nowLabel, std::chrono::steady_clock::now()) };
	writeTimeDiff(startTS, ts.get());
	timeStructs_m.push_back(std::move(ts));
}

void LogFile::writeTimeDiffFromNow(int startTSind, std::string nowLabel)
{
	writeTimeDiffFromNow(timeStructs_m.at(startTSind).get(), nowLabel);
}


void LogFile::printTimeNowFromLastTS()
{
	std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timeStructs_m.at(timeStructs_m.size()-1)->tp).count()) / 1000000;
}

void LogFile::printTimeNowFromFirstTS()
{
	std::cout << static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timeStructs_m.at(0)->tp).count()) / 1000000;
}