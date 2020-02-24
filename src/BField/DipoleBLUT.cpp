#include "BField/DipoleB.h"
#include "BField/DipoleBLUT.h"

#include <iostream>
#include <filesystem>
#include "utils/serializationHelpers.h"

using std::cerr;
using std::string;
using std::invalid_argument;
using namespace utils::fileIO::serialize;


vector<double> DipoleBLUT::getAllAttributes() const
{
	vector<double> ret{ ILAT_m, ds_msmt_m, ds_gradB_m, simMin_m, simMax_m, static_cast<double>(numMsmts_m) };
	
	for (int iii = 0; iii < altitude_m.size(); iii++)
		ret.push_back(altitude_m.at(iii));
	for (int iii = 0; iii < magnitude_m.size(); iii++)
		ret.push_back(magnitude_m.at(iii));

	return ret;
}

void DipoleBLUT::serialize(string serialFolder) const
{
	string filename{ serialFolder + string("BModel_DipoleBLUT.ser") };

	if (std::filesystem::exists(filename))
		cerr << __func__ << ": Warning: filename exists: " << filename << " You are overwriting an existing file.";
	
	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument(__func__ + string(": unable to create file: ") + filename);

	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	// ======== write data to file ======== //
	out.write(reinterpret_cast<const char*>(this), sizeof(DipoleBLUT));
	
	writeStrBuf(serializeDoubleVector(altitude_m));
	writeStrBuf(serializeDoubleVector(magnitude_m));
	
	out.close();
}

void DipoleBLUT::deserialize(string serialFolder)
{
	string filename{ serialFolder + string("/BModel_DipoleBLUT.ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument(__func__ + string(": unable to open file: ") + filename);

	DipoleBLUT* dipb{ nullptr };
	vector<char> dipbchar(sizeof(DipoleBLUT), '\0');

	in.read(dipbchar.data(), sizeof(DipoleBLUT));
	dipb = reinterpret_cast<DipoleBLUT*>(dipbchar.data());

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