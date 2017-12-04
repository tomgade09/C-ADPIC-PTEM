#include "include\LogFile.h"

void LogFile::writeLogFileEntry(std::string logData, std::string logMessage)
{//[Time (ms from start) - 15 chars tot, 11 chars for numbers] | Log Data - 15 chars | Log Message - unlimited chars
	std::string writeTxt;
	writeTxt = "[ " + std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - timeStructs_m[0].tp).count(); //time
	writeTxt += " ]"; //not sure why I can't just put this on the end of the above line but doesn't work
	
	if ((15 - writeTxt.length()) > 0)
	{
		for (int iii = 0; iii < (15 - writeTxt.length()); iii++)
			writeTxt += ' ';
	}

	writeTxt += " | " + logData;

	if ((33 - writeTxt.length()) > 0)
	{
		for (int iii = 0; iii < (33 - writeTxt.length()); iii++)
			writeTxt += ' ';
	}

	writeTxt += " | " + logMessage + "\n";

	fileIO::writeTxtFile(logFileName_m.c_str(), writeTxt.c_str());
}

void LogFile::createTimeStruct(std::string label)
{
	timeStruct tS;
	tS.label = label;
	tS.tp = std::chrono::steady_clock::now();
	timeStructs_m.push_back(tS);
}

void LogFile::writeLogTimeDiff(timeStruct startTS, timeStruct endTS)
{
	std::string logData;
	std::string logMessage;

	logData = std::chrono::duration_cast<std::chrono::milliseconds>(endTS.tp - startTS.tp).count();
	logMessage = "Time measurement " + startTS.label + " to " + endTS.label;

	writeLogFileEntry(logData, logMessage);
}
//overloads of above
void LogFile::writeLogTimeDiff(int startTSind, timeStruct endTS)
{
	writeLogTimeDiff(timeStructs_m[startTSind], endTS);
}

void LogFile::writeLogTimeDiff(int startTSind, int endTSind)
{
	writeLogTimeDiff(timeStructs_m[startTSind], timeStructs_m[endTSind] );
}


void LogFile::writeLogTimeDiffFromNow(timeStruct startTS, std::string nowLabel)
{
	timeStruct ts;
	ts.label = nowLabel;
	ts.tp = std::chrono::steady_clock::now();

	writeLogTimeDiff(startTS, ts);
}
//overloads of above
void LogFile::writeLogTimeDiffFromNow(int startTSind, std::string nowLabel)
{
	writeLogTimeDiffFromNow(timeStructs_m[startTSind], nowLabel);
}








void LogFile::printTimeNowFromTimeStruct(timeStruct* tS, std::string label)
{
	std::cout << "Time measurement " << tS->label << " to " << label << ": ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tS->tp).count() << " ms\n";
}

void LogFile::printTimeNowFromTSJustMS(timeStruct* tS)
{
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - tS->tp).count();
}

void LogFile::printTimeDiffBtwTwoTimeStructs(timeStruct* startTS, timeStruct* endTS)
{
	std::cout << "Time measurement " << startTS->label << " to " << endTS->label << ": ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTS->tp - startTS->tp).count() << " ms\n";
}

void LogFile::printTimeDiffJustTimeMS(timeStruct* startTS, timeStruct* endTS)
{
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTS->tp - startTS->tp).count();
}