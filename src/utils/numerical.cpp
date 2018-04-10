#include "utils\numerical.h"

namespace utils
{
	namespace numerical
	{
		DLLEXP_NOEXTC void normalize(std::vector<double>& normalizeMe, double normFactor, bool inverse) //inverse defaults to false
		{
			if (normFactor == 1.0)
				return;

			for (int elems = 0; elems < normalizeMe.size(); elems++) //normalize -> divide by normalization factor
				normalizeMe.at(elems) *= (inverse ? (normFactor) : (1 / normFactor));
		}

		DLLEXP_NOEXTC double calcMean(const std::vector<double>& calcMyMean, bool absValue)
		{
			double sum{ 0 };
			for (int iii = 0; iii < calcMyMean.size(); iii++)
			{
				if (absValue)
					sum += abs(calcMyMean.at(iii));
				else
					sum += calcMyMean.at(iii);
			}
			return sum / calcMyMean.size();
		}

		DLLEXP_NOEXTC double calcStdDev(std::vector<double>& calcMyStdDev)
		{
			double stdDev{ 0 };
			double mean{ calcMean(calcMyStdDev, false) };
			for (int iii = 0; iii < calcMyStdDev.size(); iii++)
			{
				stdDev += pow(calcMyStdDev.at(iii) - mean, 2);
			}
			stdDev = sqrt(stdDev / calcMyStdDev.size());
			return stdDev;
		}

		DLLEXP_NOEXTC void coutMinMaxErr(const std::vector<double>& basevals, const std::vector<double>& testvals, std::string label, bool skipzeroes) //label defaults to "", skipzeroes to true
		{
			if (basevals.size() != testvals.size())
				throw std::invalid_argument("coutMinMaxErr: vectors are not the same size");

			double maxerr{ 0.0 };
			double minerr{ 1.0e300 };
			for (int iii = 0; iii < basevals.size(); iii++)
			{
				if (basevals.at(iii) == 0.0 && skipzeroes) { continue; }
				if (testvals.at(iii) == 0.0 && skipzeroes) { continue; }
				double err{ abs((basevals.at(iii) - testvals.at(iii)) / basevals.at(iii)) };
				if (err > maxerr) { maxerr = err; }
				if (err < minerr) { minerr = err; }
			}

			std::cout << label << " min err: " << minerr << ", max err: " << maxerr << std::endl;
		}

		DLLEXP_NOEXTC void coutNumAboveErrEps(const std::vector<double>& basevals, const std::vector<double>& testvals, double errEps, std::string label, bool skipzeroes) //label defaults to "", skipzeroes to true
		{
			if (basevals.size() != testvals.size())
				throw std::invalid_argument("coutNumAboveEps: vectors are not the same size");

			int above{ 0 };
			for (int iii = 0; iii < basevals.size(); iii++)
			{
				if (basevals.at(iii) == 0.0 && skipzeroes) { continue; }
				if (testvals.at(iii) == 0.0 && skipzeroes) { continue; }
				if (abs((basevals.at(iii) - testvals.at(iii)) / basevals.at(iii)) > errEps) { above++; };
			}

			std::cout << label << " error above " << errEps << ": " << above << std::endl;
		}
	}
}