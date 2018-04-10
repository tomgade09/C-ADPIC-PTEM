#include "utils\write.h"

namespace utils
{
	namespace write
	{
		void CSV::write() //private
		{
			if (data_m.size() == 0)
			{
				printf("CSV::write: data_m empty, nothing to write to disk\n");
				return;
			}
			
			std::ofstream file(filename_m, std::ios::trunc);
			for (int iii = 0; iii < labels_m.size(); iii++)
				file << labels_m.at(iii) << ((iii != labels_m.size() - 1) ? "," : "\n");
			file.close();
			
			fileIO::write2DCSV(data_m, filename_m, (int)data_m.at(0).size(), (int)data_m.size(), ',', false);
			data_m.clear();
		}

		void ParticleDistribution::write() //private
		{
			if (data_m.at(0).size() == 0)
			{
				printf("ParticleDistribution::write: data_m empty, nothing to write to disk\n");
				return;
			}

			std::cout << "Writing particle distribution: " << std::endl;
			std::cout << saveFolder_m + "/" + particleName_m + "_" + attrNames_m.at(0) + ".bin" << std::endl;
			std::cout << saveFolder_m + "/" + particleName_m + "_" + attrNames_m.at(1) + ".bin" << std::endl;
			std::cout << saveFolder_m + "/" + particleName_m + "_" + attrNames_m.at(2) + ".bin" << std::endl;
			std::cout << data_m.size() << " " << data_m.at(0).size() << " " << data_m.at(1).size() << " " << data_m.at(2).size() << std::endl;

			for (int attr = 0; attr < data_m.size(); attr++)
				fileIO::writeDblBin(data_m.at(attr), saveFolder_m + "/" + particleName_m + "_" + attrNames_m.at(attr) + ".bin", numberParticles_m);

			data_m.clear();
			data_m = std::vector<std::vector<double>>(3);
		}

		void ParticleDistribution::addEnergyRange(int energyBins, double Emin, double Emax, bool logE) //logE defaults to true
		{
			energyPitch_m.at(0).resize(energyPitch_m.at(0).size() + energyBins);

			double binsize{ (Emax - Emin) / (energyBins - 1) };
			ranges_m.at(0).push_back({ Emin, Emax, binsize });
			
			for (int eng = 0; eng < energyBins; eng++)
				energyPitch_m.at(0).at(eng) = (logE) ? pow(10, eng * binsize + Emin) : (Emin + eng * binsize);
		}

		void ParticleDistribution::addPitchRange(int pitchBins, double PAmin, double PAmax, bool midBin) //midBin defaults to true
		{
			energyPitch_m.at(1).resize(energyPitch_m.at(1).size() + pitchBins);
			
			double binsize{ (midBin) ? ((PAmax - PAmin) / pitchBins) : ((PAmax - PAmin) / (pitchBins - 1)) };
			ranges_m.at(1).push_back({ PAmin, PAmax, binsize });

			for (int ang = 0; ang < pitchBins; ang++)
				energyPitch_m.at(1).at(ang) = (midBin) ? (PAmax - (ang + 0.5) * binsize) : (PAmax - ang * binsize);
		}

		void ParticleDistribution::generate(double s_ion, double s_mag)
		{
			numberParticles_m = (int)(energyPitch_m.at(0).size() * energyPitch_m.at(1).size());
			std::vector<double> s;

			for (auto pit = energyPitch_m.at(1).begin(); pit < energyPitch_m.at(1).end(); pit++)
				for (int eng = 0; eng < energyPitch_m.at(0).size(); eng++) { s.push_back((*pit) <= 90 ? s_mag : s_ion); }

			generate(s);
		}

		void ParticleDistribution::generate(std::vector<double>& s)
		{
			if (numberParticles_m == 0) { numberParticles_m = (int)(energyPitch_m.at(0).size() * energyPitch_m.at(1).size()); }

			std::cout << "Energy Bins (min, max : bin size): ";
			for (int eng = 0; eng < ranges_m.at(0).size(); eng++)
				std::cout << ((eng > 0) ? "                                   " : "") << ranges_m.at(0).at(eng).at(0) << ", " << ranges_m.at(0).at(eng).at(1) << " : " << ranges_m.at(0).at(eng).at(2) << std::endl;
			std::cout << "Pitch Bins (min, max : bin size): ";
			for (int ang = 0; ang < ranges_m.at(1).size(); ang++)
				std::cout << ((ang > 0) ? "                                  " : "") << ranges_m.at(1).at(ang).at(0) << ", " << ranges_m.at(1).at(ang).at(1) << " : " << ranges_m.at(1).at(ang).at(2) << std::endl;

			for (auto pit = energyPitch_m.at(1).begin(); pit < energyPitch_m.at(1).end(); pit++)
			{
				for (auto eng = energyPitch_m.at(0).begin(); eng < energyPitch_m.at(0).end(); eng++)
				{
					data_m.at(0).push_back(-sqrt(2 * (*eng) * JOULE_PER_EV / mass_m) * cos((*pit) * PI / 180.0));
					data_m.at(1).push_back( sqrt(2 * (*eng) * JOULE_PER_EV / mass_m) * sin((*pit) * PI / 180.0));
				}
			}

			data_m.at(2) = s;

			std::cout << "First two vpara/vperp: " << data_m.at(0).at(0) << "/" << data_m.at(1).at(0) << "  " << data_m.at(0).at(1) << "/" << data_m.at(1).at(1) << std::endl;
		}
	} /* end namespace utils::write */
} /* end namespace utils */