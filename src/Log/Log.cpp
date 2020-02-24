#include <iostream>
#include <iomanip>
#include <sstream>

#include "Log/Log.h"

using std::cout;
using std::cerr;
using std::setw;
using std::to_string;
using std::make_unique;
using std::setprecision;
using std::stringstream;
using std::invalid_argument;

using namespace std::chrono;

Log::Entry::Entry(std::chrono::steady_clock::time_point time, string message, bool write, bool error) :
	time_m{ time }, message_m{ message }, write_m{ write }, error_m{ error }
{

}

void Log::save()
{
	ofstream logFile;

	logFile.open(logFilePath_m);
	logFile << string("[  Time (ms)  ] : Log Message\n");

	auto writeEntry = [&](const Entry& entry)
	{
		if (!entry.write_m) return;

		duration<double, std::nano> dur = entry.time_m - entries_m.at(0).time_m;

		string writeTxt;
		stringstream timess;

		timess << setprecision(10) << setw(11) << dur.count() / 1000000;
		string timestr{ timess.str() };
		while (timestr.length() > 11)
			timestr.erase(std::ios::end);
		writeTxt = "[ " + timestr + " ] : " + ((entry.error_m) ? "ERROR: " : "") + entry.message_m + "\n";

		logFile << writeTxt;
	};

	for (const auto& entr : entries_m)
		writeEntry(entr);

	logFile.close();
}

void Log::cerrCheck()
{
	while (check_m)
	{
		if (cerr_m.str().length() > 0)
		{
			size_t tmplen{ cerr_m.str().length() };

			do std::this_thread::sleep_for(microseconds(100));
			while (tmplen != cerr_m.str().length()); //block until cerr is not being written to anymore

			entries_m.push_back(Entry(steady_clock::now(), cerr_m.str(), true, true));
			cerr_m.str(string());
			cerr_m.clear();
		}

		std::this_thread::sleep_for(milliseconds(10)); //check cerr every __ <<
	}
}

Log::Log(string logFilePath) : logFilePath_m{ logFilePath }
{
	cerr.rdbuf(cerr_m.rdbuf()); //set cerr read buffer to this class' stringstream

	cerrReadThread_m = std::thread([&]() { cerrCheck(); });

	createEntry("Simulation initialized");
}

Log::~Log()
{
	check_m = false;
	cerrReadThread_m.join();

	cerr.rdbuf(cerrBufferBackup_m); //restore cerr to normal
	save();
}

size_t Log::createEntry(string label, bool write)
{
	entries_m.push_back(Entry(steady_clock::now(), label, write));

	return entries_m.size() - 1;
}

double Log::timeElapsedTotal_s() const
{
	return timeElapsedSinceEntry_s(0);
}

double Log::timeElapsedSinceEntry_s(size_t index) const
{
	duration<double, std::nano> elapsed = steady_clock::now() - entries_m.at(index).time_m;
	return elapsed.count() / 1000000000.0; //convert to s
}