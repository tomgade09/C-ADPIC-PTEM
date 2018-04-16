#ifndef SIMULATIONATTRIBUTES_H
#define SIMULATIONATTRIBUTES_H

#include <vector>
#include <string>
#include <functional>
#include "FileIO\fileIO.h"
#include "utils\string.h"

#define VEC(x) std::vector<x> //to save lots of space
#define STR std::string //to save lots of space

class SimAttributes
{
private:
	struct attrsData
	{
		const std::string classname_m;
		VEC(STR)          names_m;
		VEC(VEC(STR))     strLabels_m;
		VEC(VEC(STR))     strAttrs_m;
		VEC(VEC(STR))     dblLabels_m;
		VEC(VEC(double))  dblAttrs_m;

		attrsData(std::string classname) : classname_m{ classname } {}
	};

	std::string saveString_m;
	std::string filename_m;
	bool read_m{ false };

	std::string generateString(attrsData& ad);
	void write();
	void read();

	std::function<attrsData*(std::string)> matchClassname = [&](std::string check) {
		if (check == "Simulation") return &simAD;
		else if (check == "BField") return &BAD;
		else if (check == "EField")	return &EAD;
		else if (check == "Particle") return &partAD;
		else if (check == "Satellite") return &satAD;
		else throw std::invalid_argument("SimAttributes::addData: invalid argument - no class of name " + check); };


public: //all this is public so callers can access the raw data (so I don't have to write equally many access functions which doesn't make sense)
	attrsData simAD{ "Simulation" }; //Simulation attributes
	attrsData BAD{ "BField" }; //BField attributes
	attrsData EAD{ "EField" }; //EField attributes
	attrsData partAD{ "Particle" }; //Particle attributes
	attrsData satAD{ "Satellite" }; //Satellite attributes

	SimAttributes(std::string filename, bool readFile = false) : filename_m{ filename }, read_m{ readFile } { if (readFile) { read(); read_m = true; } }
	~SimAttributes() { if (!read_m) write(); }
	
	void addData(STR classname, STR name, VEC(STR) stringAttrLabels, VEC(STR) stringAttributes, VEC(STR) doubleAttrLabels, VEC(double) doubleAttributes);
};

#undef VEC
#undef STR

#endif /* !SIMULATIONATTRIBUTES_H */