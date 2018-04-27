#include "test.h"

#include <vector>

namespace test
{
	//functions not visible outside the file
	std::vector<double> dipoleGPU()
	{
		return std::vector<double>();
	}

	std::vector<double> dipoleLUTGPU()
	{
		return std::vector<double>();
	}

	//outside accessible tests
	bool dipoleAccuracy()
	{
		return true;
	}

	bool QSPSAccuracy()
	{
		return true;
	}
}