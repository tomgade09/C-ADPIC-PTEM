#ifndef UTILS_LOAD_H
#define UTILS_LOAD_H

#include "fileIO\FileIO.h"
#include "utils\numerical.h"
#include "utils\string.h"

namespace utils
{
	namespace load
	{
		class DistributionFromDisk
		{
		private:
			std::vector<std::string> attrNames_m;
			std::vector<std::vector<double>> data_m;


		public:
			DistributionFromDisk(std::string folder, std::string partName, std::vector<std::string> attrNames) : attrNames_m{ attrNames }
			{
				int attrsize{ 0 };
				for (int attr = 0; attr < attrNames.size(); attr++)
				{
					std::vector<double> read;
					fileIO::readDblBin(read, folder + "/" + partName + "_" + attrNames.at(attr) + ".bin");
					data_m.push_back(read);
					if (attrNames_m.at(attr).size() > attrsize) { attrsize = attrNames_m.at(attr).size(); }
				}

				for (int attr = 0; attr < attrNames.size(); attr++)
					if (attrNames_m.at(attr).size() < attrsize) { utils::string::stringPadder(attrNames_m.at(attr), attrsize); }
			}
			~DistributionFromDisk() {}

			const std::vector<std::vector<double>>& data() { return data_m; }
			void print(int at) {
				for (int iii = 0; iii < attrNames_m.size(); iii++) { std::cout << attrNames_m.at(iii) << ((iii != attrNames_m.size() - 1) ? ", " : ":"); }
				for (int iii = 0; iii < data_m.size(); iii++) { std::cout << data_m.at(iii).at(at) << ((iii != data_m.size() - 1) ? ", " : ""); } std::cout << std::endl; }
			void printdiff(DistributionFromDisk& other, int at);
			void zeroes() { std::vector<int> tmp; zeroes(tmp, true); }
			void zeroes(std::vector<int>& zeroes, bool print = true);
			void compare(DistributionFromDisk& other);
		};
	}
}

#endif /* UTILS_LOAD_H */