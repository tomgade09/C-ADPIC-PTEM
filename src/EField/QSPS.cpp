#include "EField/QSPS.h"

#include <iostream>
#include <filesystem>

#include "utils/serializationHelpers.h"

using std::cerr;
using std::to_string;
using std::invalid_argument;
using namespace utils::fileIO::serialize;

void QSPS::serialize(string serialFolder) const
{
	string filename{ serialFolder + string("EField_QSPS.ser") };

	if (std::filesystem::exists(filename))
		cerr << "QSPS::serialize: Warning: filename exists: " << filename << " You are overwriting an existing file.";

	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument("QSPS::serialize: unable to create file: " + filename);

	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	// ======== write data to file ======== //
	out.write(reinterpret_cast<const char*>(this), sizeof(QSPS));
	writeStrBuf(serializeString(string(name_m)));
	writeStrBuf(serializeDoubleVector(altMin_m));
	writeStrBuf(serializeDoubleVector(altMax_m));
	writeStrBuf(serializeDoubleVector(magnitude_m));

	out.close();
}

void QSPS::deserialize(string serialFolder, int nameIndex)
{
	string filename{ serialFolder + string("EField_QSPS" + to_string(nameIndex) + ".ser") };

	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument("QSPS::deserialize: unable to open file: " + filename);

	QSPS* qsps{ nullptr };
	vector<char> qspschar(sizeof(QSPS));

	in.read(qspschar.data(), sizeof(QSPS));
	qsps = reinterpret_cast<QSPS*>(qspschar.data());

	name_m = deserializeString(in).c_str();
	altMin_m = deserializeDoubleVector(in);
	altMax_m = deserializeDoubleVector(in);
	magnitude_m = deserializeDoubleVector(in);

	useGPU_m = qsps->useGPU_m;

	this_d = nullptr;
}