#include "utils/writeIOclasses.h"
#include "utils/fileIO.h"

#include <iostream>
#include <fstream>
#include <cmath>

namespace utils
{
	namespace fileIO
	{
		//start CSV private functions
		void CSV::write() //private
		{
			if (data_m.size() == 0)
			{
				printf("CSV::write: data_m empty, nothing to write to disk\n");
				return;
			}

			std::ofstream file(filename_m, std::ios::trunc);
			for (unsigned int iii = 0; iii < labels_m.size(); iii++)
				file << labels_m.at(iii) << ((iii != labels_m.size() - 1) ? "," : "\n");
			file.close();

			fileIO::write2DCSV(data_m, filename_m, (int)data_m.at(0).size(), (int)data_m.size(), ',', false);
			data_m.clear();
		}

		void CSV::pad(std::vector<double>& vec)
		{
			if (vec.size() > data_m.size())
				printf("CSV::pad: warning: vec is bigger than data_m.at(0) - some doubles in vec will be trimmed off\n");
			vec.resize(data_m.size());
		}

		//public
		CSV::~CSV()
		{
			write();
		}

		void CSV::add(std::vector<double> vec, std::string label)
		{
			if (data_m.size() > 0 && vec.size() != data_m.at(0).size())
			{
				pad(vec);
			}
			data_m.push_back(vec);
			labels_m.push_back(label);
		}

		void CSV::add(std::vector<std::vector<double>> vecs, std::vector<std::string> labels)
		{
			if (vecs.size() != labels.size())
				throw std::invalid_argument("CSV::add: vecs.size() != labels.size()");
			for (auto vec = vecs.begin(); vec < vecs.end(); vec++)
				add((*vec), labels.at(vec - vecs.begin()));
		}

		void CSV::addspace()
		{
			if (data_m.size() == 0)
			{
				printf("CSV::addspace: data_m is empty, cannot add zeroes array\n");
				return;
			}
			data_m.push_back(std::vector<double>(data_m.at(0).size()));
			labels_m.push_back("");
		}


		//start ParticleDistribution private member functions
		double ParticleDistribution::EPAtoV(double E_eV, double PA_deg, bool vpara)
		{
			if (E_eV > 0.0)
			{
				if (vpara)
					return -sqrt(2 * E_eV * JOULE_PER_EV / mass_m) * cos(PA_deg * RADS_PER_DEG);
				else
					return  sqrt(2 * E_eV * JOULE_PER_EV / mass_m) * sin(PA_deg * RADS_PER_DEG);
			}
			else if (E_eV == 0.0)
				return 0.0;
			else
				throw std::invalid_argument("ParticleDistribution::EPAtoV: negative energy used as argument - this is invalid/unphysical - E(eV): " + std::to_string(E_eV));
		}

		void ParticleDistribution::write()
		{
			if (data_m.at(0).size() == 0)
			{
				printf("ParticleDistribution::write: data_m empty, nothing to write to disk\n");
				return;
			}

			std::cout << "Writing particle distribution: " << std::endl;
			for (auto& attrNm : attrNames_m)
				std::cout << saveFolder_m + "/" + particleName_m + "_" + attrNm + ".bin" << std::endl;
			std::cout << "[" << data_m.size() << "] x [" << data_m.at(0).size() << "]" << std::endl;

			for (size_t attr = 0; attr < data_m.size(); attr++)
				fileIO::writeDblBin(data_m.at(attr), saveFolder_m + "/" + particleName_m + "_" + attrNames_m.at(attr) + ".bin", (unsigned int)data_m.at(0).size());

			clear();
			data_m = std::vector<std::vector<double>>(3);
		}

		//public
		ParticleDistribution::ParticleDistribution(std::string saveFolder, std::vector<std::string> attrNames, std::string particleName, double mass, std::vector<double> padvals) :
			saveFolder_m{ saveFolder }, attrNames_m{ attrNames }, particleName_m{ particleName }, mass_m{ mass }, padvals_m{ padvals }
		{
			ranges_m = std::vector<std::vector<std::vector<double>>>(2);
			energyPitch_m = std::vector<std::vector<double>>(2);
			data_m = std::vector<std::vector<double>>(attrNames_m.size());
		}

		ParticleDistribution::~ParticleDistribution()
		{
			write();
		}

		void ParticleDistribution::setattr(std::vector<double>& attr, unsigned int ind)
		{
			if (ind >= attrNames_m.size())
				throw std::invalid_argument("ParticleDistribution::setattr: ind is higher than data_m.size() - 1");
			data_m.at(ind) = attr;
		}

		void ParticleDistribution::setattr(std::vector<double>& attr, std::string name)
		{
			for (unsigned int name_ind = 0; name_ind < attrNames_m.size(); name_ind++)
				if (name == attrNames_m.at(name_ind)) { setattr(attr, name_ind); return; }

			throw std::logic_error("ParticleDistribution::setattr: no attribute named " + name);
		}

		void ParticleDistribution::addEnergyRange(unsigned int energyBins, double E_start, double E_end, bool logE) //logE defaults to true
		{
			int bins{ (int)energyBins };
			int oldsize{ (int)energyPitch_m.at(0).size() };
			energyPitch_m.at(0).resize(oldsize + bins);

			double binsize{ (E_end - E_start) / (bins - 1) };
			ranges_m.at(0).push_back({ E_start, E_end, binsize });

			for (int eng = 0; eng < bins; eng++) //starts at E_start and increments to E_end
				energyPitch_m.at(0).at(oldsize + eng) = (logE) ? pow(10, E_start + eng * binsize) : (E_start + eng * binsize);
		}

		void ParticleDistribution::addPitchRange(unsigned int pitchBins, double PA_start, double PA_end, bool midBin) //midBin defaults to true
		{
			int bins{ (int)pitchBins };
			int oldsize{ (int)energyPitch_m.at(1).size() };
			energyPitch_m.at(1).resize(oldsize + bins);

			double binsize{ (midBin) ? ((PA_end - PA_start) / bins) : ((PA_end - PA_start) / (bins - 1)) };
			ranges_m.at(1).push_back({ PA_start, PA_end, binsize });

			for (int ang = 0; ang < bins; ang++) //starts at PA_start and increments to PA_end
				energyPitch_m.at(1).at(oldsize + ang) = (midBin) ? (PA_start + (ang + 0.5) * binsize) : (PA_start + ang * binsize);
		}

		void ParticleDistribution::addSpecificParticle(unsigned int numParticles, double energy, double pitch, double s, unsigned int padmult)
		{
			if (ranges_m.at(0).size() != 0 || ranges_m.at(1).size() != 0)
			{
				std::cout << "ParticleDistribution::addSpecificParticle: At least one of { energy, pitch } has at least one range specified.  addSpecificParticle is not (for now) compatible with specifying ranges.  Clearing existing ranges." << std::endl;
				clear(); //this is because data_m is generated by iterating each energy at each pitch - handling/not modifying existing data in the array would require some workarounds I haven't coded yet.
			}

			double vpara{ EPAtoV(energy, pitch, true ) };
			double vperp{ EPAtoV(energy, pitch, false) };

			unsigned int finalsize{ (padmult != 0) ?
				((((unsigned int)data_m.at(0).size() + numParticles) / padmult) + 1) * padmult : //calculate the nearest larger binsize that is a multiple of padmult
				(unsigned int)data_m.at(0).size() + numParticles //if not padding, just calculate current size + number to add
			};
			unsigned int padnum{ finalsize - ((unsigned int)data_m.at(0).size() + numParticles) }; //how many placeholders do we have to add to get to the above size? - 0 if not padding, something else if padding

			for (auto& attr : data_m)
				attr.reserve(finalsize);

			for (unsigned int part = 0; part < numParticles; part++)
			{
				data_m.at(0).push_back(vpara);
				data_m.at(1).push_back(vperp);
				data_m.at(2).push_back(s);

				for (unsigned int otherattr = 3; otherattr < data_m.size(); otherattr++)
					data_m.at(otherattr).push_back(padvals_m.at(otherattr)); //pad the remaining attributes
			}

			if (padmult != 0)
			{//used to pad the number of values so length is even multiple of padmult (for CUDA execution - make a multiple of blocksize)
				for (unsigned int attr = 0; attr < data_m.size(); attr++)
					for (unsigned int pad = 0; pad < padnum; pad++)
						data_m.at(attr).push_back(padvals_m.at(attr));
			}
		}

		void ParticleDistribution::generate(double s_ion, double s_mag)
		{
			std::vector<double> s;

			for (auto& pitch : energyPitch_m.at(1))
				for (auto& energy : energyPitch_m.at(0))
					s.push_back(pitch <= 90.0 ? s_mag : s_ion);

			generate(s);
		}

		void ParticleDistribution::generate(std::vector<double>& s)
		{
			if (energyPitch_m.at(0).size() == 0 || energyPitch_m.at(1).size() == 0)
			{
				std::cout << "ParticleDistribution::addSpecificParticle: At least one of { energy, pitch } range is not specified.  Returning." << std::endl;
				return;
			}

			std::cout << "Energy Bins (min, max : bin size): ";
			for (unsigned int eng = 0; eng < ranges_m.at(0).size(); eng++)
				std::cout << ((eng > 0) ? "                                   " : "") << ranges_m.at(0).at(eng).at(0) << ", " << ranges_m.at(0).at(eng).at(1) << " : " << ranges_m.at(0).at(eng).at(2) << std::endl;
			std::cout << "Pitch Bins (min, max : bin size): ";
			for (unsigned int ang = 0; ang < ranges_m.at(1).size(); ang++)
				std::cout << ((ang > 0) ? "                                  " : "") << ranges_m.at(1).at(ang).at(0) << ", " << ranges_m.at(1).at(ang).at(1) << " : " << ranges_m.at(1).at(ang).at(2) << std::endl;

			unsigned int finalsize{ ((unsigned int)energyPitch_m.at(0).size() * (unsigned int)energyPitch_m.at(1).size() + (unsigned int)data_m.at(0).size()) }; //energy bins * pitch bins + any other particles currently specified = total

			for (auto& attr : data_m)
				attr.reserve(finalsize); //prevents multiple resize/deep copy operations

			for (auto& pitch : energyPitch_m.at(1))
			{
				for (auto& energy : energyPitch_m.at(0))
				{
					data_m.at(0).push_back(EPAtoV(energy, pitch, true ));
					data_m.at(1).push_back(EPAtoV(energy, pitch, false));

					for (unsigned int otherattr = 3; otherattr < data_m.size(); otherattr++)
						data_m.at(otherattr).push_back(padvals_m.at(otherattr)); //pad the remaining attributes
				}
			}

			data_m.at(2) = s;
		}

		void ParticleDistribution::clear()
		{
			ranges_m.clear();
			data_m.clear();
			energyPitch_m.clear();

			ranges_m = std::vector<std::vector<std::vector<double>>>(2);
			energyPitch_m = std::vector<std::vector<double>>(2);
			data_m = std::vector<std::vector<double>>(attrNames_m.size());
		}
	} /* end namespace utils::write */
} /* end namespace utils */
