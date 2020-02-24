#include "EField/EField.h"
#include "EField/QSPS.h"

#include <iostream>
#include <sstream>
#include <filesystem>

#include "utils/serializationHelpers.h"

using std::cerr;
using std::move;
using std::string;
using std::to_string;
using std::make_unique;
using std::stringstream;
using std::runtime_error;
using std::invalid_argument;
using namespace utils::fileIO::serialize;

void EField::serialize(string serialFolder) const
{
	string filename{ serialFolder + string("/EField.ser") };

	if (std::filesystem::exists(filename))
		cerr << "EField::serialize: Warning: filename exists: " << filename << " You are overwriting an existing file.\n";

	ofstream out(filename, std::ofstream::binary);
	if (!out) throw invalid_argument("EField::serialize: unable to create file: " + filename);

	auto writeStrBuf = [&](const stringbuf& sb)
	{
		out.write(sb.str().c_str(), sb.str().length());
	};

	out.write(reinterpret_cast<const char*>(this), sizeof(EField));
	
	size_t numModels{ emodels_m.size() };
	out.write(reinterpret_cast<char*>(&numModels), sizeof(size_t));

	for (const auto& emodel : emodels_m)
	{
		out.write(reinterpret_cast<char*>(&(emodel->type_m)), sizeof(EModel::Type)); //write type of emodel
		writeStrBuf(emodel->serialize());
	}

	out.close();
}

void EField::deserialize(string serialFolder)
{
	string filename{ serialFolder + string("/EField.ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument("EField::deserialize: unable to open file: " + filename);

	vector<char> efieldchar(sizeof(EField), '\0');
	in.read(reinterpret_cast<char*>(efieldchar.data()), sizeof(EField));
	
	useGPU_m = (*reinterpret_cast<EField*>(efieldchar.data())).useGPU_m;
	
	size_t len{ readSizetLength(in) };
	
	for (size_t emodel = 0; emodel < len; emodel++)
	{
		EModel::Type type{ -1 };
		in.read(reinterpret_cast<char*>(&type), sizeof(EModel::Type));

		if (type == EModel::Type::QSPS) emodels_m.push_back(move(make_unique<QSPS>(in)));
		//else if (type == EModel::Type::AlfvenLUT) elements_m.push_back(move(make_unique<AlfvenLUT>(in)));
		else throw runtime_error("EField::deserialize: unknown EModel Type");
	}
}