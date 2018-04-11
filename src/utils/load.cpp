#include "utils\load.h"

namespace utils
{
	namespace load
	{
		void DistributionFromDisk::printdiff(DistributionFromDisk& other, int at)
		{
			int datasize{ (int)data_m.size() };
			if (data_m.size() != other.data().size())
			{
				std::cout << "DistributionFromDisk::printdiff: Warning: data from the two distributions does not have the same dimensionality. ";
				std::cout << "Did you load up two distributions of different types? Using the smaller size." << std::endl;
				datasize = (data_m.size() < other.data().size()) ? ((int)data_m.size()) : ((int)other.data().size());
			}

			std::vector<double> err(datasize);
			for (int attr = 0; attr < datasize; attr++)
			{
				err.at(attr) = abs((data_m.at(attr).at(at) - other.data().at(attr).at(at)) / data_m.at(attr).at(at));
				std::cout << attrNames_m.at(attr) << " (this, other, err): " << data_m.at(attr).at(at) << ", " << other.data().at(attr).at(at) << ", " << err.at(attr) << std::endl;
			}
		}

		void DistributionFromDisk::zeroes(std::vector<int>& zeroes, bool print) //print defaults to true
		{
			if (zeroes.size() != data_m.size()) { zeroes.resize(data_m.size()); }
			for (auto attr = data_m.begin(); attr < data_m.end(); attr++)
			{
				for (auto part = (*attr).begin(); part < (*attr).end(); part++)
					if ((*part) == 0.0) { zeroes.at(attr - data_m.begin())++; }
			}

			for (int zero = 0; zero < zeroes.size(); zero++)
				std::cout << attrNames_m.at(zero) << " zeroes: " << zeroes.at(zero) << (zero != zeroes.size() - 1) ? ", " : "";
			std::cout << std::endl;
		}

		void DistributionFromDisk::compare(DistributionFromDisk& other)
		{
			const std::vector<std::vector<double>>& data_other{ other.data() };
			int datasize{ (int)data_m.size() };
			if (data_m.size() != data_other.size())
			{
				std::cout << "DistributionFromDisk::compare: Warning: data from the two distributions does not have the same dimensionality. ";
				std::cout << "Did you load up two distributions of different types? Using the smaller size." << std::endl;
				datasize = (data_m.size() < data_other.size()) ? ((int)data_m.size()) : ((int)data_other.size());
			}

			std::vector<int> zeroes_this;
			std::vector<int> zeroes_other;
			zeroes(zeroes_this, false);
			other.zeroes(zeroes_other, false);

			std::vector<double> mean_this;
			std::vector<double> mean_other;
			std::vector<double> stdDev_this;
			std::vector<double> stdDev_other;

			for (int iii = 0; iii < datasize; iii++)
			{
				mean_this.push_back(utils::numerical::calcMean(data_m.at(iii)));
				mean_other.push_back(utils::numerical::calcMean(data_other.at(iii)));
				stdDev_this.push_back(utils::numerical::calcStdDev(data_m.at(iii)));
				stdDev_other.push_back(utils::numerical::calcStdDev(data_other.at(iii)));
			}

			std::vector<int> notsame(datasize);
			std::vector<double> avgErr(datasize);
			std::vector<double> minErr(datasize);
			std::vector<double> maxErr(datasize);
			for (int attr = 0; attr < datasize; attr++)
			{
				int partsize{ (int)data_m.at(attr).size() };
				if (data_m.at(attr).size() != data_other.at(attr).size()) { std::cout << "DistributionFromDisk::compare: Warning: attributes have different number of particles.  Using smaller number" << std::endl;
					partsize = (data_m.at(attr).size() < data_other.at(attr).size()) ? ((int)data_m.at(attr).size()) : ((int)data_other.at(attr).size());	}
				for (int part = 0; part < partsize; part++)
				{
					if (data_m.at(attr).at(part) != data_other.at(attr).at(part))
					{
						double err{ abs((data_m.at(attr).at(part) - data_other.at(attr).at(part)) / data_m.at(attr).at(part)) };
						if (minErr.at(attr) > err) { minErr.at(attr) = err; }
						if (maxErr.at(attr) < err) { maxErr.at(attr) = err; }
						avgErr.at(attr) = (avgErr.at(attr) * notsame.at(attr) + err) / (notsame.at(attr) + 1);
						notsame.at(attr)++;
					}
				}
			}

			std::cout << "================ Summary of differences: this, other (passed in) ================" << std::endl;
			std::cout << "Number of zeroes in attributes:" << std::endl;
			for (int attr = 0; attr < datasize; attr++)
				std::cout << attrNames_m.at(attr) << ": " << zeroes_this.at(attr) << ", " << zeroes_other.at(attr) << std::endl;
			std::cout << std::endl;

			std::cout << "Attribute means:" << std::endl;
			for (int attr = 0; attr < datasize; attr++)
				std::cout << attrNames_m.at(attr) << ": " << mean_this.at(attr) << ", " << mean_other.at(attr) << std::endl;
			std::cout << std::endl;

			std::cout << "Attribute standard deviations:" << std::endl;
			for (int attr = 0; attr < datasize; attr++)
				std::cout << attrNames_m.at(attr) << ": " << stdDev_this.at(attr) << ", " << stdDev_other.at(attr) << std::endl;
			std::cout << std::endl;

			std::cout << "Attribute error (other's deviation from this):" << std::endl;
			for (int attr = 0; attr < datasize; attr++)
				std::cout << attrNames_m.at(attr) << ": Min: " << minErr.at(attr) << ", Max: " << maxErr.at(attr) << ", Avg: " << avgErr.at(attr) << ", Number not same: " << notsame.at(attr) << std::endl;
			std::cout << std::endl;
		}
	}
}