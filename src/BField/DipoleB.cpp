#include "BField/DipoleB.h"

#include <iostream>
#include <filesystem>
#include "utils/serializationHelpers.h"

using std::cerr;
using std::string;
using std::invalid_argument;
using namespace utils::fileIO::serialize;

vector<double> DipoleB::getAllAttributes() const
{
	return { L_m,L_norm_m, s_max_m, ILAT_m, ds_m, lambdaErrorTolerance_m };
}

void DipoleB::serialize(string serialFolder) const
{
	string filename{ serialFolder + string("BModel_DipoleB.ser") };

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

	out.close();
}

void DipoleB::deserialize(string serialFolder)
{
	string filename{ serialFolder + string("/BModel_DipoleB.ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument(__func__ + string(": unable to open file: ") + filename);

	DipoleB* dipb{ nullptr };
	vector<char> dipbchar(sizeof(DipoleB), '\0');

	in.read(dipbchar.data(), sizeof(DipoleB));
	dipb = reinterpret_cast<DipoleB*>(dipbchar.data());

	this_d = nullptr;

	L_m = dipb->L_m;
	L_norm_m = dipb->L_norm_m;
	s_max_m = dipb->s_max_m;
	ILAT_m = dipb->ILAT_m;
	ds_m = dipb->ds_m;
	lambdaErrorTolerance_m = dipb->lambdaErrorTolerance_m;

	useGPU_m = dipb->useGPU_m;
}