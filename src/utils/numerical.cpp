#include "utils/numerical.h"
#include <cmath>

namespace utils
{
	namespace numerical
	{
		DLLEXP void v2DtoEPitch(const std::vector<double>& vpara, const std::vector<double>& vperp, double mass, std::vector<double>& energies_eV, std::vector<double>& pitches_deg)
		{
			if (vpara.size() != vperp.size())
				throw std::invalid_argument("utils::numerical::v2DtoEPitch: input vectors vpara and vperp are not the same size: " + std::to_string(vpara.size()) + ", " + std::to_string(vperp.size()));

			if (energies_eV.size() != vpara.size()) //resize output vectors to be as big as input
				energies_eV.resize(vpara.size());
			if (pitches_deg.size() != vpara.size())
				pitches_deg.resize(vpara.size());

			std::fill(energies_eV.begin(), energies_eV.end(), 0.0); //zeroes data in the return arrays
			std::fill(pitches_deg.begin(), pitches_deg.end(), 0.0); //guarantees E/PA is 0 if vpara/vperp is zero

			for (size_t part = 0; part < vpara.size(); part++)
			{
				bool nonZero{ vpara.at(part) != 0.0 || vperp.at(part) != 0.0 };

				if (nonZero) //check this or else the function can produce "NaN" in some indicies (I think atan2 is responsible) -> if false, the data at that index will be left 0
				{
					energies_eV.at(part) = (0.5 * mass * (vpara.at(part) * vpara.at(part) + vperp.at(part) * vperp.at(part)) / JOULE_PER_EV);
					pitches_deg.at(part) = atan2(std::abs(vperp.at(part)), -vpara.at(part)) / RADS_PER_DEG;
				}
			}
		}

		DLLEXP void EPitchTov2D(const std::vector<double>& energies_eV, const std::vector<double>& pitches_deg, double mass, std::vector<double>& vpara, std::vector<double>& vperp)
		{
			if (energies_eV.size() != pitches_deg.size())
				throw std::invalid_argument("utils::numerical::EPitchTov2D: input vectors vpara and vperp are not the same size: " + std::to_string(vpara.size()) + ", " + std::to_string(vperp.size()));

			if (energies_eV.size() != vpara.size()) //resize output vectors to be as big as input
				vpara.resize(energies_eV.size());
			if (energies_eV.size() != vperp.size())
				vperp.resize(energies_eV.size());

			std::fill(vpara.begin(), vpara.end(), 0.0); //zeroes data in the return arrays
			std::fill(vperp.begin(), vperp.end(), 0.0); //guarantees E/PA is 0 if vpara/vperp is zero

			for (size_t part = 0; part < energies_eV.size(); part++)
			{
				bool nonZero{ energies_eV.at(part) != 0.0 };

				if (nonZero) //check this or else the function can produce "NaN" in some indicies (I think atan2 is responsible) -> if false, the data at that index will be left 0
				{
					vpara.at(part) = -sqrt(2 * energies_eV.at(part) * JOULE_PER_EV / mass) * cos(pitches_deg.at(part) * RADS_PER_DEG);
					vperp.at(part) =  sqrt(2 * energies_eV.at(part) * JOULE_PER_EV / mass) * sin(pitches_deg.at(part) * RADS_PER_DEG);
				}
			}
		}

		DLLEXP std::vector<double> generateSpacedValues(double start, double end, int number, bool logSpaced, bool endInclusive)
		{
			/*
				**Note** if logSpaced is true, min and max have to be log(min) and log(max),
				min/max and in between values will be used as powers of 10 (10^(min | max))

				x -> values that are returned
			     dval
			    |-----|
				-------------------------
				x     x     x     x     x
				-------------------------
				^start ============= end^
				^min                 max^ << endInclusive = true, ("number" - 1) values
				^min           max^       << endInclusive = false, "number" values
			*/

			if (number <= 0)
				throw std::invalid_argument("utils::numerical::generateSpacedValues: number of values is less than / equal to zero");

			std::vector<double> ret(number);

			double dval{ (end - start) / ((endInclusive) ? (number - 1) : number) };
			for (int iter = 0; iter < number; iter++)
				ret.at(iter) = ((logSpaced) ? pow(10, iter * dval + start) : (iter * dval + start));

			return ret;
		}

		DLLEXP void normalize(std::vector<double>& normalizeMe, double normFactor, bool inverse) //inverse defaults to false
		{
			if (normFactor == 1.0)
				return;

			for (auto& elem : normalizeMe) //normalize -> divide by normalization factor
				elem *= (inverse ? (normFactor) : (1 / normFactor));
		}

		DLLEXP double calcMean(const std::vector<double>& calcMyMean, bool absValue) //absValue defaults to false
		{
			double sum{ 0 };
			for (size_t iii = 0; iii < calcMyMean.size(); iii++)
			{
				if (absValue)
					sum += std::abs(calcMyMean.at(iii));
				else
					sum += calcMyMean.at(iii);
			}
			return sum / calcMyMean.size();
		}

		DLLEXP double calcStdDev(const std::vector<double>& calcMyStdDev)
		{
			double stdDev{ 0 };
			double mean{ calcMean(calcMyStdDev, false) };
			for (size_t iii = 0; iii < calcMyStdDev.size(); iii++)
			{
				stdDev += pow(calcMyStdDev.at(iii) - mean, 2);
			}
			stdDev = sqrt(stdDev / calcMyStdDev.size());
			return stdDev;
		}

		DLLEXP void coutMinMaxErr(const std::vector<double>& basevals, const std::vector<double>& testvals, std::string label, bool skipzeroes) //label defaults to "", skipzeroes to true
		{
			if (basevals.size() != testvals.size())
				throw std::invalid_argument("coutMinMaxErr: vectors are not the same size");

			double maxerr{ 0.0 };
			double minerr{ 1.0e300 };
			for (size_t iii = 0; iii < basevals.size(); iii++)
			{
				if (basevals.at(iii) == 0.0 && skipzeroes) { continue; }
				if (testvals.at(iii) == 0.0 && skipzeroes) { continue; }
				double err{ std::abs((basevals.at(iii) - testvals.at(iii)) / basevals.at(iii)) };
				if (err > maxerr) { maxerr = err; }
				if (err < minerr) { minerr = err; }
			}

			std::cout << label << " min err: " << minerr << ", max err: " << maxerr << std::endl;
		}

		DLLEXP void coutNumAboveErrEps(const std::vector<double>& basevals, const std::vector<double>& testvals, double errEps, std::string label, bool skipzeroes) //label defaults to "", skipzeroes to true
		{
			if (basevals.size() != testvals.size())
				throw std::invalid_argument("coutNumAboveEps: vectors are not the same size");

			int above{ 0 };
			for (size_t iii = 0; iii < basevals.size(); iii++)
			{
				if (basevals.at(iii) == 0.0 && skipzeroes) { continue; }
				if (testvals.at(iii) == 0.0 && skipzeroes) { continue; }
				if (std::abs((basevals.at(iii) - testvals.at(iii)) / basevals.at(iii)) > errEps) { above++; };
			}

			std::cout << label << " error above " << errEps << ": " << above << std::endl;
		}
	}
}
