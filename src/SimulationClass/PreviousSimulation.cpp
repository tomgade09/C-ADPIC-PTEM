#include "SimulationClass\PreviousSimulation.h"

int PreviousSimulation::findAttrInd(std::string attr, std::vector<std::string> allAttrs)
{
	for (int ind = 0; ind < allAttrs.size(); ind++)
	{
		if (allAttrs.at(ind) == attr)
			return ind;
	}

	std::string allAttrsStr;
	for (int attr = 0; attr < allAttrs.size(); attr++)
		allAttrsStr += allAttrs.at(attr);

	throw std::invalid_argument("PreviousSimulation::findAttrInd: cannot find attribute " + attr + " in string " + allAttrsStr);
}