#include "SimAttributes\SimAttributes.h"

#include <iostream>

#include "utils\fileIO.h"
#include "utils\string.h"

using utils::fileIO::writeTxtFile;
using utils::fileIO::readTxtFile;
using utils::string::strToStrVec;

/*
	Saving paradigm:
	[[[[ ]]]] encapsulates classnames
	(((( )))) encapsulates name_m's from classes
	{{{{ }}}} encapsulate data following
	<<<< >>>> encapsulates all data for class
	;;;; separates entries
	$$$$ #$^& encapsulates double values converted to string

	example: <<<<[[[[Particle]]]]((((elec)))){{attrLbl,attrLbl,attrLbl,attrLbl}}{{vpara,vperp,s,tincd}}{{numparticles,somethingelse}}{{$$$$doublehere$$$$,$$$$doublehere$$$$}};;;;
		((elecbs)){{strhere,strhere}}{{strattr,strattr}}{{strdbl,strdbl}}{{$$$$doublehere$$$$,$$$$doublehere$$$$}}>>>>
*/

#define VEC(x) std::vector<x> //to save lots of space
void SimAttributes::addData(std::string classname, std::string name, VEC(std::string) stringAttrLabels, VEC(std::string) stringAttributes, VEC(std::string) doubleAttrLabels, VEC(double) doubleAttributes)
{
	auto invalidCharChk = [](std::vector<std::string> vec) { for (auto str = vec.begin(); str < vec.end(); str++)
		{ for (auto chr = (*str).begin(); chr < (*str).end(); chr++) { std::string tmp{ (*chr) };
			if (tmp == "[" || tmp == "]" || tmp == "(" || tmp == ")" || tmp == "{" || tmp == "}" ||
				tmp == "<" || tmp == ">" || tmp == "$" || tmp == "#" || tmp == ";")
				throw std::invalid_argument("SimAttributes::addData: invalid character in a string vector: " + tmp); } } };

	//function guards
	if ((stringAttrLabels.size() != stringAttributes.size()) || (doubleAttrLabels.size() != doubleAttributes.size()))
		{ throw std::invalid_argument("SimAttributes::addData: an attribute vector and its labels are not identical in size"); }
	invalidCharChk({ name });
	invalidCharChk(stringAttrLabels);
	invalidCharChk(stringAttributes);
	invalidCharChk(doubleAttrLabels);

	if (classname == "Simulation" && simAD.names_m.size() != 0) { throw std::invalid_argument("SimAttributes::addData: Simulation has already been specified"); }
	if (classname == "BField" && BAD.names_m.size() != 0) { throw std::invalid_argument("SimAttributes::addData: BField has already been specified"); }

	attrsData* adptr{ matchClassname(classname) };

	adptr->names_m.push_back(name);
	adptr->strLabels_m.push_back(stringAttrLabels);
	adptr->strAttrs_m.push_back(stringAttributes);
	adptr->dblLabels_m.push_back(doubleAttrLabels);
	adptr->dblAttrs_m.push_back(doubleAttributes);
}

std::string SimAttributes::generateString(attrsData& ad)
{
	auto dblToExactStr = [](double d) { std::string dblstr; dblstr.resize(8);
		for (int iii = 0; iii < 8; iii++) { dblstr[iii] = reinterpret_cast<char*>(&d)[iii]; } return dblstr; };

	auto strVec1DToStr = [](VEC(std::string)& strvec) {
		std::string strout;
		for (auto attr = strvec.begin(); attr < strvec.end(); attr++)
			strout += (*attr) + ((attr != strvec.end() - 1) ? "," : "");
		return strout; };

	#define ALLDATWRAP(x) std::string("<<<<") + x + std::string(">>>>")
	#define CLASSWRAP(x) std::string("[[[[") + x + std::string("]]]]")
	#define NAMEWRAP(x) std::string("((((") + x + std::string("))))")
	#define ATTRSWRAP(x) std::string("{{{{") + x + std::string("}}}}")
	#define DBLSTRWRAP(x) std::string("$$$$") + x + std::string("#$^&")
	std::string ret;

	ret = CLASSWRAP(ad.classname_m);
	for (int entry = 0; entry < ad.strLabels_m.size(); entry++) //iterates over entries
	{
		ret += NAMEWRAP(ad.names_m.at(entry));
		ret += ATTRSWRAP(strVec1DToStr(ad.strLabels_m.at(entry)));
		ret += ATTRSWRAP(strVec1DToStr(ad.strAttrs_m.at(entry)));
		ret += ATTRSWRAP(strVec1DToStr(ad.dblLabels_m.at(entry)));
		
		std::string tmp;
		for (auto dbl = ad.dblAttrs_m.at(entry).begin(); dbl < ad.dblAttrs_m.at(entry).end(); dbl++)
			tmp += DBLSTRWRAP(dblToExactStr(*dbl)) + ((dbl != ad.dblAttrs_m.at(entry).end() - 1) ? "," : "");

		ret += ATTRSWRAP(tmp);
		if (entry != ad.strLabels_m.size() - 1) { ret += ";;;;"; }
	}

	return ALLDATWRAP(ret);
}

void SimAttributes::write()
{
	saveString_m = generateString(simAD); //clear string in case anything is there
	saveString_m += generateString(BAD);
	saveString_m += generateString(EAD);
	saveString_m += generateString(partAD);
	saveString_m += generateString(satAD);

	writeTxtFile(saveString_m, filename_m, true);
}

void SimAttributes::read()
{
	readTxtFile(saveString_m, filename_m);
	
	//Liberal use of lambda functions
	auto findCutString = [](std::string between, std::string& findstr, bool erase = true) {
		std::string front{ between.substr(0,between.size() / 2) }; std::string back{ between.substr(between.size() / 2, between.size() / 2) };
		auto beg{ findstr.find(front) + front.size() }; auto end{ findstr.find(back) + back.size() };
		if (end - beg == 0) { return std::string(""); }
		std::string ret{ findstr.substr(beg, end - between.size()) }; if (erase) { findstr.erase(0, end); } return ret; };

	auto splitEntries = [](std::string splitchars, std::string& findstr) { std::vector<std::string> ret;
		while (findstr.find(splitchars) != std::string::npos)
			{ ret.push_back(findstr.substr(0, findstr.find(splitchars))); findstr.erase(0, findstr.find(splitchars) + 4); }
		ret.push_back(findstr); findstr.clear();
		return ret; };

	auto entrToAtrVec = [&](std::vector<std::string>& entryStrings) { std::vector<std::vector<std::vector<std::string>>> ret(5);
		for (auto entr = entryStrings.begin(); entr < entryStrings.end(); entr++) {
			ret.at(0).push_back( { findCutString(NAMEWRAP(""), (*entr)) } ); //name
			ret.at(1).push_back( strToStrVec(findCutString(ATTRSWRAP(""), (*entr)), ',') ); //strLabels
			ret.at(2).push_back( strToStrVec(findCutString(ATTRSWRAP(""), (*entr)), ',') ); //strAttrs
			ret.at(3).push_back( strToStrVec(findCutString(ATTRSWRAP(""), (*entr)), ',') ); //dblLabels
			ret.at(4).push_back( { findCutString(ATTRSWRAP(""), (*entr)) } ); //double string - don't want to search for ',' in case that's in the double string
			
			std::cout << (*entr) << std::endl;
		}
		return ret; };

	auto exactStr2Dbl = [](std::string s) { double ret; for (int iii = 0; iii < 8; iii++) { reinterpret_cast<char*>(&ret)[iii] = s[iii]; } return ret; };

	auto strvToDblVec = [&](std::vector<std::vector<std::string>>& vec) { std::vector<std::vector<double>> ret;
		for (int iii = 0; iii < vec.size(); iii++)
		{
			std::vector<double> tmp;
			while (vec.at(iii).at(0).find("$$$$") != std::string::npos)
				tmp.push_back(exactStr2Dbl(findCutString(DBLSTRWRAP(""), vec.at(iii).at(0))));
			ret.push_back(tmp);
		}
		return ret;	};
	//end lambdas - just for this function! (finally)

	std::string saveString{ saveString_m };

	std::vector<std::string> dataStrVec; //vector that holds data Classes at indicies - Simulation, BField...etc
	dataStrVec.push_back(findCutString(ALLDATWRAP(""), saveString)); //Simulation
	dataStrVec.push_back(findCutString(ALLDATWRAP(""), saveString)); //BField
	dataStrVec.push_back(findCutString(ALLDATWRAP(""), saveString)); //EField
	dataStrVec.push_back(findCutString(ALLDATWRAP(""), saveString)); //Particle
	dataStrVec.push_back(findCutString(ALLDATWRAP(""), saveString)); //Satellite

	for (auto datastr = dataStrVec.begin(); datastr < dataStrVec.end(); datastr++) //datastr is a data string that contains all data for a class
	{ //iterates and pushes each Class to various functions resulting in attrsData having all the loaded data
		attrsData* adptr{ matchClassname(findCutString(CLASSWRAP(""), (*datastr))) };
		std::vector<std::vector<std::vector<std::string>>> allDataStr{ entrToAtrVec(splitEntries(";;;;", (*datastr))) };

		for (auto name = allDataStr.at(0).begin(); name < allDataStr.at(0).end(); name++)
			adptr->names_m.push_back((*name).at(0));
		adptr->strLabels_m = allDataStr.at(1);
		adptr->strAttrs_m = allDataStr.at(2);
		adptr->dblLabels_m = allDataStr.at(3);
		adptr->dblAttrs_m = strvToDblVec(allDataStr.at(4));
	}

	if (EAD.names_m.size() == 1 && EAD.names_m.at(0) == "")
	{
		EAD.names_m.clear();
		EAD.strLabels_m.clear();
		EAD.strAttrs_m.clear();
		EAD.dblLabels_m.clear();
		EAD.dblAttrs_m.clear();
	}
}

#undef VEC
#undef ALLDATWRAP
#undef CLASSWRAP
#undef NAMEWRAP
#undef ATTRSWRAP
#undef DBLSTRWRAP