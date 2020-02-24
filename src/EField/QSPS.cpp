#include "EField/QSPS.h"

#include <iostream>
#include <filesystem>

using std::cerr;
using std::to_string;
using std::invalid_argument;
using namespace utils::fileIO::serialize;

vector<double> QSPS::getAllAttributes() const
{
	vector<double> ret;
	for (int iii = 0; iii < altMin_m.size(); iii++)
	{//vectors are guaranteed to be the same size
		ret.push_back(altMin_m.at(iii));
		ret.push_back(altMax_m.at(iii));
		ret.push_back(magnitude_m.at(iii));
	}

	return ret;
}

stringbuf QSPS::serialize() const
{
	stringbuf sb;
	ostream out(&sb);

	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	// ======== write data to file ======== //
	//out.write(reinterpret_cast<char*>(type_m), sizeof(Type)); //written by EField
	out.write(reinterpret_cast<const char*>(this), sizeof(QSPS));
	writeStrBuf(serializeDoubleVector(altMin_m));
	writeStrBuf(serializeDoubleVector(altMax_m));
	writeStrBuf(serializeDoubleVector(magnitude_m));

	return sb;
}

void QSPS::deserialize(ifstream& in)
{
	vector<char> typechar(sizeof(Type), '\0');
	in.read(typechar.data(), sizeof(Type));

	vector<char> qspschar(sizeof(QSPS), '\0');
	in.read(qspschar.data(), sizeof(QSPS));

	altMin_m = deserializeDoubleVector(in);
	altMax_m = deserializeDoubleVector(in);
	magnitude_m = deserializeDoubleVector(in);

	useGPU_m = (*reinterpret_cast<QSPS*>(qspschar.data())).useGPU_m;
	
	this_d = nullptr;
}