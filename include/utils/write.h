#ifndef UTILS_WRITE_H
#define UTILS_WRITE_H

#include <vector>
#include <string>

#include "physicalconstants.h"

namespace utils
{
	namespace write
	{
		class CSV
		{
		private:
			std::string filename_m;
			
			std::vector<std::vector<double>> data_m;
			std::vector<std::string> labels_m;

			void write(); //defined in cpp

			void pad(std::vector<double>& vec) { if (vec.size() > data_m.size()) { printf("CSV::pad: warning: vec is bigger than data_m.at(0) - some doubles in vec will be trimmed off\n"); }
				vec.resize(data_m.size()); }

		public:
			CSV(std::string filename) : filename_m{ filename } {}
			~CSV() { write(); }

			void add(std::vector<double> vec, std::string label) { if (data_m.size() > 0 && vec.size() != data_m.at(0).size())
				{ pad(vec); } data_m.push_back(vec); labels_m.push_back(label); }

			void add(std::vector<std::vector<double>> vecs, std::vector<std::string> labels) { if (vecs.size() != labels.size()) { throw std::invalid_argument("CSV::add: vecs.size() != labels.size()"); }
				for (auto vec = vecs.begin(); vec < vecs.end(); vec++) { add((*vec), labels.at(vec - vecs.begin())); } }
			
			void addspace() { if (data_m.size() == 0) { printf("CSV::addspace: data_m is empty, cannot add zeroes array\n"); return; }
				data_m.push_back(std::vector<double>(data_m.at(0).size())); labels_m.push_back(""); }

			std::vector<std::vector<double>>& data() { return data_m; }
		};

		class ParticleDistribution
		{
		private:
			std::string saveFolder_m;
			std::string particleName_m;
			std::vector<std::string> attrNames_m;

			double mass_m;

			std::vector<std::vector<std::vector<double>>> ranges_m; //[energy || pitch][range index][min, max, step]
			std::vector<std::vector<double>> energyPitch_m; //[energy || pitch][energy/pitch index]
			std::vector<std::vector<double>> data_m; //[vpara || vperp || s (|| time || index)][particle index]

			double EPAtoV(double E_eV, double PA_deg, bool vpara);
			void write();


		public:
			ParticleDistribution(std::string saveFolder = "./", std::vector<std::string> attrNames = { "vpara", "vperp", "s" }, std::string particleName = "elec", double mass = MASS_ELECTRON):
				saveFolder_m{ saveFolder }, attrNames_m{ attrNames }, particleName_m{ particleName }, mass_m{ mass }
			{
				ranges_m = std::vector<std::vector<std::vector<double>>>(2);
				energyPitch_m = std::vector<std::vector<double>>(2);
				data_m = std::vector<std::vector<double>>(attrNames_m.size());
			}
			~ParticleDistribution() { write(); }

			const std::vector<std::vector<double>>& data() { return data_m; }
			void setattr(std::vector<double>& attr, int ind) { if (ind >= attrNames_m.size()) { throw std::invalid_argument("ParticleDistribution::setattr: ind is higher than data_m.size() - 1"); }
				data_m.at(ind) = attr; }

			void addEnergyRange(int energyBins, double Emin, double Emax, bool logE = true);
			void addPitchRange(int pitchBins, double PAmin, double PAmax, bool midBin = true);
			void addSpecificParticle(int numParticles, double energy, double pitch, double s, int padmult = 0);
			void generate(double s_ion, double s_mag);
			void generate(std::vector<double>& s);

			void clear();
		};
	}
}

#endif /* !UTILS_WRITE_H */