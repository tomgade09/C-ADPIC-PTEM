#ifndef UTILS_IOCLASSES_READ_H
#define UTILS_IOCLASSES_READ_H

#include <vector>
#include <string>

#include "dlldefines.h"

namespace utils
{
	namespace fileIO
	{
		class DistributionFromDisk
		{
		private:
			std::string name_m;
			std::vector<std::string> attrNames_m;
			std::vector<std::vector<double>> data_m;


		public:
			DistributionFromDisk(std::string name, std::string folder, std::string partName, std::vector<std::string> attrNames);
			~DistributionFromDisk() {}

			const std::vector<std::vector<double>>& data() const { return data_m; }
			const std::string& name() const { return name_m; }
			void print(int at) const;
			void printdiff(DistributionFromDisk& other, int at) const;
			void zeroes() const;
			void zeroes(std::vector<int>& zeroes, bool print = true) const;
			void compare(const DistributionFromDisk& other) const;
		};
	}
}

#endif /* UTILS_IOCLASSES_READ_H */