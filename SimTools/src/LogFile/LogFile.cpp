#include "LogFile\LogFile.h"

void LogFile::writeLogFileEntry(std::string logData, std::string logMessage)
{//[Time (ms from start) - 19 chars tot, 15 chars for numbers] | Log Data - 20 chars | Log Message - unlimited chars
	std::string writeTxt;
	writeTxt = "[ " + std::to_string(static_cast<float>(std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - timeStructs_m[0]->tp).count())/1000) + " ]"; //time

	size_t txtlen = writeTxt.length();

	if ((19 - txtlen) > 0)
	{
		for (int iii = 0; iii < (19 - txtlen); iii++)
			writeTxt += ' ';
	}
	else
		writeTxt.erase(17, txtlen - 19);

	writeTxt += " | " + logData;
	txtlen = writeTxt.length();

	if ((42 - txtlen) > 0)
	{
		for (int iii = 0; iii < (42 - txtlen); iii++)
			writeTxt += ' ';
	}
	else
		writeTxt.erase(42, txtlen - 42);

	writeTxt += " | " + logMessage + "\n";
	fileIO::writeTxtFile(logFileName_m.c_str(), writeTxt.c_str());
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
	std::string logData;
	std::string logMessage;

	logData = "writeTimeDiff";
	logMessage = "Timed:  " + startTS->label + "  TO  " + endTS->label + ": " \
		+ std::to_string(static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(endTS->tp - startTS->tp).count()) / 1000);

	writeLogFileEntry(logData, logMessage);
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