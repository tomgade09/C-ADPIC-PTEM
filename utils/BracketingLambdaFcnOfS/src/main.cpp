#include <cmath>
#include <iostream>
#include <iomanip>

#define DLLFILE
#ifdef DLLFILE
#define DLLEXPORT extern "C" __declspec(dllexport)
#else
#define DLLEXPORT
#endif

constexpr double PI{ 3.14159265359 };
constexpr double ERRTOLERANCE{ 1e-4 };
constexpr double B0{ 3.12e-5 }; //in T, the Earth's magnetic field at the equator and surface of the Earth, used for the dipole equations
constexpr double Re{ 6.371e6 };
bool firstrun{ true };

DLLEXPORT double getSAtLambda(double lambdaDegrees, double L)///FIX TO GIVE S FROM RE NOT EQUATOR!!!!!!!!!!!AA!!!1111!1!!111!
{//returns s in units of L
	double x{ std::asinh(std::sqrt(3.0) * std::sin(lambdaDegrees * PI / 180.0)) };

	return (0.5 * L / std::sqrt(3.0)) * (x + std::sinh(x) * std::cosh(x));
}

DLLEXPORT double getLambdaAtS(double s, double dipoleEquatorialDist, double ILATdegrees, double dlambda0, double lambDivFact, int* numIters)
{//s, L, and dipoleEquatorialDist must be in same units
	double s_max{ getSAtLambda(ILATdegrees, dipoleEquatorialDist) };
	double lambda_tmp{ (-ILATdegrees / s_max) * s + ILATdegrees };
	double s_tmp{ s_max - getSAtLambda(lambda_tmp, dipoleEquatorialDist) };
	double dlambda{ dlambda0 };
	bool   over{ 0 };
	int    iters{ 0 };

	//if (dlambda0 == 1.5 && lambDivFact == 1.25) { std::cout << std::endl << lambda_tmp << "  " << s_tmp << "  " << dlambda << "  " << lambDivFact << std::endl; }

	while (std::abs((s_tmp - s) / s) > ERRTOLERANCE)
	{
		while (1)
		{
			iters++;
			over = (s_tmp >= s);
			if (over)
			{
				lambda_tmp += dlambda;
				s_tmp = s_max - getSAtLambda(lambda_tmp, dipoleEquatorialDist);
				if (s_tmp < s)
					break;
			}
			else
			{
				lambda_tmp -= dlambda;
				s_tmp = s_max - getSAtLambda(lambda_tmp, dipoleEquatorialDist);
				if (s_tmp >= s)
					break;
			}
		}
		dlambda /= lambDivFact;
	}

	*numIters = iters;

	return lambda_tmp;
}

double BFieldAtS(double L, double s_max, double ILAT, double s)
{// consts: [ B0, ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	int tmp{ 0 };
	double lambda_deg{ getLambdaAtS(s, L, ILAT, 1.0, 5.0, &tmp) };
	double rnorm{ L / Re * cos(PI * lambda_deg / 180.0) * cos(PI * lambda_deg / 180.0) };

	return -3.12e-5 / (rnorm * rnorm * rnorm) * sqrt(1.0 + 3 * sin(PI * lambda_deg / 180.0) * sin(PI * lambda_deg / 180.0));
}

double gradBAtS(double L, double s_max, double ILAT, double s, double ds = 6371.2)
{
	return (BFieldAtS(L, s_max, ILAT, s + ds) - BFieldAtS(L, s_max, ILAT, s - ds)) / (2 * ds);
}

int main()
{
	double sIn{ 0.0 };
	double ILAT{ 72.0 };
	double L{ Re / pow(cos(ILAT * PI / 180.0), 2) };
	double s_max{ getSAtLambda(ILAT, L) };
	double ds{ 6371.2 };

	//std::cout << "Specify s in meters: ";
	//std::cin >> sIn;

	double sBtm{ 2030837.49610366 }; //bottom
	double sTop{ 19881647.2473464 }; //top

	sIn = sBtm;

	std::cout << std::setprecision(10) << sIn << "  " << ILAT << "  " << L << "  " << s_max << "  " << ds << std::endl;
	std::cout << std::setprecision(10) << "BFieldAtS: " << BFieldAtS(L, s_max, ILAT, sIn) << std::endl;

	while (sIn < sTop)
	{
		int minIters{ 100000 };
		int maxIters{ 0 };
		double best_dlambda{ 0.0 };
		double wrst_dlambda{ 0.0 };
		double best_lmbdivfact{ 0.0 };
		double wrst_lmbdivfact{ 0.0 };

		double dlambda{ 1.5 };
		double lambDivFact{ 1.25 };

		for (int iii = 0; iii < 15; iii++)
		{
			for (int jjj = 0; jjj < 36; jjj++)
			{
				//dlambda = 1.5 - 0.1 * iii;
				//lambDivFact = 1.25 + 0.25 * jjj;

				dlambda = 1.2;
				lambDivFact = 3.0;

				int numIters{ 0 };
				getLambdaAtS(sIn, L, ILAT, dlambda, lambDivFact, &numIters);

				if (minIters > numIters)
				{
					minIters = numIters;
					best_dlambda = dlambda;
					best_lmbdivfact = lambDivFact;
				}
				if (maxIters < numIters)
				{
					maxIters = numIters;
					wrst_dlambda = dlambda;
					wrst_lmbdivfact = lambDivFact;
				}
			}
		}

		std::cout << std::setprecision(6) << "s (Re)," << sIn / Re << ",Lowest iters:," << minIters << ", dlmb:," << best_dlambda << ", lmb/x:," << best_lmbdivfact << ",|| Most iters:," << maxIters << ", dlmb:," << wrst_dlambda << ", lmb/x:," << wrst_lmbdivfact << std::endl;

		sIn += Re / 10;
	}
	
	/*double minusds{ 500 }; //ds for grad B accuracy test
	double plusds{ 50000 };

	for (int iii = 0; iii < 13; iii++)
	{
		std::cout << std::setprecision(10) << "gradBAtS," << ds - iii * minusds << "," << gradBAtS(L, s_max, ILAT, sIn, ds - iii * minusds) << std::endl;
	}
	std::cout << "plus" << std::endl;
	for (int iii = 0; iii < 50; iii++)
	{
		std::cout << std::setprecision(10) << "gradBAtS," << ds + iii * plusds << "," << gradBAtS(L, s_max, ILAT, sIn, ds - iii * plusds) << std::endl;
	}*/

	return 0;
}