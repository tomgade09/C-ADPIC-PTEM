#include "BField/DipoleB.h"
#include "BField/DipoleBLUT.h"

#include <iostream>
#include <filesystem>
#include "utils/serializationHelpers.h"

using std::cerr;
using std::string;
using std::invalid_argument;
using namespace utils::fileIO::serialize;


void DipoleBLUT::serialize(string serialFolder) const
{
	string filename{ serialFolder + string("BField_DipoleBLUT.ser") };

	if (std::filesystem::exists(filename))
		cerr << __func__ << ": Warning: filename exists: " << filename << " You are overwriting an existing file.";
	
	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument(__func__ + string(": unable to create file: ") + filename);

	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	// ======== write data to file ======== //
	out.write(reinterpret_cast<const char*>(this), sizeof(DipoleB));
	writeStrBuf(serializeString(string(name_m)));
	
	serializeDoubleVector(altitude_m);
	serializeDoubleVector(magnitude_m);
	
	out.close();
}

void DipoleBLUT::deserialize(string serialFolder)
{
	string filename{ serialFolder + string("/BField_DipoleBLUT.ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument(__func__ + string(": unable to open file: ") + filename);

	DipoleBLUT* dipb{ nullptr };
	vector<char> dipbchar(sizeof(DipoleBLUT));

	in.read(dipbchar.data(), sizeof(DipoleBLUT));
	dipb = reinterpret_cast<DipoleBLUT*>(dipbchar.data());

	name_m = deserializeString(in).c_str();
	this_d = nullptr;

	ILAT_m = dipb->ILAT_m;
	//ds_msmt_m = dipb->ds_msmt_m; //is calculated in the ctor
	ds_gradB_m = dipb->ds_gradB_m;

	altitude_m = deserializeDoubleVector(in);
	magnitude_m = deserializeDoubleVector(in);
	altitude_d = nullptr;
	magnitude_d = nullptr;

	simMin_m = dipb->simMin_m;
	simMax_m = dipb->simMax_m;
	numMsmts_m = dipb->numMsmts_m;

	useGPU_m = dipb->useGPU_m;
}