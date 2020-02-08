#include "EField/EField.h"
#include "EField/QSPS.h"

#include <iostream>
#include <sstream>
#include <filesystem>

#include "utils/serializationHelpers.h"

using std::cerr;
using std::string;
using std::to_string;
using std::stringstream;
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

	size_t size{ Eelems_m.size() };
	out.write(reinterpret_cast<char*>(&size), sizeof(size_t));

	int QSPScnt{ 0 };
	int ALUTcnt{ 0 };
	for (const auto& elem : Eelems_m)
	{
		if (elem->name() == "QSPS") ++QSPScnt;//writeStrBuf(serializeString(string(elem.name()) + to_string(++QSPScnt)));
		else if (elem->name() == "AlfvenLUT") ++ALUTcnt;//writeStrBuf(serializeString(string(elem.name()) + to_string(++ALUTcnt)));
		else throw invalid_argument("EField::serialize: element does not have a recognized name: " + string(elem->name()));
	}

	out.write(reinterpret_cast<char*>(&QSPScnt), sizeof(int));
	out.write(reinterpret_cast<char*>(&ALUTcnt), sizeof(int));

	for (const auto& elem : Eelems_m)
	{
		elem->serialize(serialFolder);
		//need to verify use of filename::exists and filename::move below
		if (elem->name() == "QSPS" && std::filesystem::exists("EField_QSPS.ser"))
		{
			int iter{ 0 };
			while (std::filesystem::exists("EField_QSPS" + to_string(iter) + ".ser")) iter++;
			std::filesystem::rename("EField_QSPS.ser", "EField_QSPS" + to_string(iter) + ".ser");
		}
		else if (elem->name() == "AlfvenLUT" && std::filesystem::exists("EField_AlfvenLUT.ser"))
		{
			int iter{ 0 };
			while (std::filesystem::exists("EField_AlfvenLUT" + to_string(iter) + ".ser")) iter++;
			std::filesystem::rename("EField_AlfvenLUT.ser", "EField_AlfvenLUT" + to_string(iter) + ".ser");
		}
	}

	out.close();
}

void EField::deserialize(string serialFolder)
{
	string filename{ serialFolder + string("/EField.ser") };
	ifstream in(filename, std::ifstream::binary);
	if (!in) throw invalid_argument("EField::deserialize: unable to open file: " + filename);

	size_t elemcnt{ readSizetLength(in) };

	int QSPScnt{ 0 };
	int ALUTcnt{ 0 };

	in.read(reinterpret_cast<char*>(&QSPScnt), sizeof(int));
	in.read(reinterpret_cast<char*>(&ALUTcnt), sizeof(int));

	for (size_t elem = 0; elem < QSPScnt; elem++)
		Eelems_m.push_back(std::move(std::make_unique<QSPS>(serialFolder, elem)));

	//for (size_t elem = 0; elem < ALUTcnt; elem++)
		//Eelems_m.push_back(std::move(std::make_unique<AlfvenLUT>(serialFolder, elem)));
}