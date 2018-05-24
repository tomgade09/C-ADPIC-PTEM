#ifndef SIMFATALEXCEPTION_H
#define SIMFATALEXCEPTION_H

#include "ErrorHandling/SimException.h"

//for signalling exit()

class SimFatalException : public SimException
{
protected:

public:
	SimFatalException(std::string error, std::string fileName, int lineNum, std::vector<double> numArgs = {}, std::vector<std::string> strArgs = {}) :
		SimException(error, fileName, lineNum, numArgs, strArgs) {}
};

#endif /* SIMFATALEXCEPTION */