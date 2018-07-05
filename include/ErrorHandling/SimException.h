#ifndef SIMULATIONERROR_H
#define SIMULATIONERROR_H

#include <exception>
#include <string>
#include <vector>

class SimException : public std::exception
{
protected:
	std::string error_m;
	std::string fileLine_m;
	std::string args_m;
	std::vector<double> numericalArguments_m;
	std::vector<std::string> stringArguments_m;

public:
	SimException(std::string error, std::string fileName, int lineNum, std::vector<double> numArgs = {}, std::vector<std::string> strArgs = {}) :
		error_m{ error }, fileLine_m{ fileName + ":" + std::to_string(lineNum) }, numericalArguments_m { numArgs }, stringArguments_m{ strArgs }
	{
		args_m = "args::: ";
		if (!numericalArguments_m.empty())
		{
			args_m += "double: ";
			for (int iii = 0; iii < numericalArguments_m.size(); iii++)
			{
				args_m += std::to_string(numericalArguments_m.at(iii));
				if (iii < numericalArguments_m.size() - 1)
					args_m += ", ";
			}
		}

		if (!stringArguments_m.empty())
		{
			args_m += "; string: ";
			for (int iii = 0; iii < stringArguments_m.size(); iii++)
			{
				args_m += stringArguments_m.at(iii);
				if (iii < stringArguments_m.size() - 1)
					args_m += ", ";
			}
		}
	}

	//~SimException() throw() {}

	virtual const char* what()  const throw() { return error_m.c_str(); }
	virtual const char* where() const { return fileLine_m.c_str(); }
	virtual const char* args()  const { return args_m.c_str(); }
};


#endif /* SIMULATIONERROR_H */
