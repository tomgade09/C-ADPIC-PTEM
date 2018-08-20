#ifndef UTILS_IOCLASSES_WRITE_H
#define UTILS_IOCLASSES_WRITE_H

#include <vector>
#include <string>
#include <stdexcept>

#include "dlldefines.h"
#include "physicalconstants.h"

namespace utils
{
	namespace fileIO
	{
		class CSV
		{
		private:
			std::string filename_m;
			
			std::vector<std::vector<double>> data_m;
			std::vector<std::string> labels_m;

			void write(); //defined in cpp
			void pad(std::vector<double>& vec);

		public:
			CSV(std::string filename) : filename_m{ filename } {}
			~CSV();

			void add(std::vector<double> vec, std::string label);
			void add(std::vector<std::vector<double>> vecs, std::vector<std::string> labels);			
			void addspace();

			std::vector<std::vector<double>>& data() { return data_m; }
		};

		class ParticleDistribution
		{
		private:
			std::string saveFolder_m;
			std::vector<std::string> attrNames_m;
			std::string particleName_m;

			double mass_m;

			std::vector<std::vector<std::vector<double>>> ranges_m; //[energy || pitch][range index][min, max, step]
			std::vector<std::vector<double>> energyPitch_m; //[energy || pitch][energy/pitch index]
			std::vector<double> padvals_m;
			std::vector<std::vector<double>> data_m; //[vpara || vperp || s (|| time || index)][particle index]

			double EPAtoV(double E_eV, double PA_deg, bool vpara);
			void write();


		public: //generate is dependent on vpara, vperp, and s being the first three attributes - if not, this will have to be modified
			ParticleDistribution(std::string saveFolder = "./", std::vector<std::string> attrNames = { "vpara", "vperp", "s", "t_inc", "t_esc" }, std::string particleName = "elec", double mass = MASS_ELECTRON, std::vector<double> padvals = { 0.0, 0.0, 0.0, 0.0, -1.0});
			~ParticleDistribution(); //writes on destruction
			ParticleDistribution(const ParticleDistribution&) = delete;
			ParticleDistribution& operator=(const ParticleDistribution&) = delete;

			const std::vector<std::vector<double>>& data() const { return data_m; }
			const std::vector<std::vector<std::vector<double>>>& ranges() const { if (ranges_m.at(0).size() == 0 || ranges_m.at(1).size() == 0) throw std::logic_error("ParticleDistribution::ranges(): one or both ranges (PA, E) not specified"); return ranges_m; }
			void setattr(std::vector<double>& attr, unsigned int ind);
			void setattr(std::vector<double>& attr, std::string name);
			void addEnergyRange(unsigned int energyBins, double E_start, double E_end, bool logE = true);
			void addPitchRange(unsigned int pitchBins, double PA_start, double PA_end, bool midBin = true);
			void addSpecificParticle(unsigned int numParticles, double energy, double pitch, double s, unsigned int padmult = 0);
			void generate(double s_ion, double s_mag);
			void generate(std::vector<double>& s);
			void clear();
		};
	}
}

#endif /* !UTILS_IOCLASSES_WRITE_H */