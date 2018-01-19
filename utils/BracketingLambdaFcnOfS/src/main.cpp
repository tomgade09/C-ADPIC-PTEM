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
constexpr double ERRTOLERANCE{ 1e-10 };
constexpr double B0{ 3.12e-5 }; //in T, the Earth's magnetic field at the equator and surface of the Earth, used for the dipole equations
constexpr double Re{ 6.371e6 };
bool firstrun{ true };

DLLEXPORT double getSAtLambda(double lambdaDegrees, double L)///FIX TO GIVE S FROM RE NOT EQUATOR!!!!!!!!!!!AA!!!1111!1!!111!
{//returns s in units of L
	double x{ std::asinh(std::sqrt(3) * std::sin(lambdaDegrees * PI / 180)) };

	return (0.5 * L / std::sqrt(3)) * (x + std::sinh(x) * std::cosh(x));
}

DLLEXPORT double getLambdaAtS(double s, double dipoleEquatorialDist, double ILATdegrees)
{//s, L, and dipoleEquatorialDist must be in same units
	double s_max{ getSAtLambda(ILATdegrees, dipoleEquatorialDist) };
	double lambda_tmp{ (-ILATdegrees / s_max) * s + ILATdegrees };
	double s_tmp{ s_max - getSAtLambda(lambda_tmp, dipoleEquatorialDist) };
	double dlambda{ 1 };
	bool   over{ 0 };

	//std::cout << "s: " << s << ", L: " << dipoleEquatorialDist << std::endl;
	//std::cout << "s_max: " << s_max / 6.371e6 << " lambda_tmp: " << lambda_tmp << " s_tmp: " << s_tmp / 6.371e6 << std::endl;

	while (std::abs((s_tmp - s) / s) > ERRTOLERANCE)
	{
		while (1)
		{
			over = (s_tmp >= s);
			if (over)
			{
				lambda_tmp += dlambda;
				s_tmp = s_max - getSAtLambda(lambda_tmp, dipoleEquatorialDist);
				//std::cout << std::setprecision(10) << "over: " << s_tmp << ", " << lambda_tmp << std::endl;
				if (s_tmp < s)
					break;
			}
			else
			{
				lambda_tmp -= dlambda;
				s_tmp = s_max - getSAtLambda(lambda_tmp, dipoleEquatorialDist);
				//std::cout << std::setprecision(10) << "under: " << s_tmp << ", " << lambda_tmp << std::endl;
				if (s_tmp >= s)
					break;
			}
		}
		dlambda /= 10;
		//std::cout << std::setprecision(10) << "dlambda: " << dlambda << std::endl;
	}

	return lambda_tmp;
}

DLLEXPORT double getBx(double s, double dipoleEquatorialDist, double ILATdegrees)
{
	if (firstrun) { std::cout << "getBx: " << s << "  " << dipoleEquatorialDist << "  " << ILATdegrees << std::endl; }
	double L{ Re / std::pow(std::cos(ILATdegrees * PI / 180), 2) };
	double lambda_final{ getLambdaAtS(s, L, ILATdegrees) };
	double lambda_final_rad{ lambda_final * PI / 180 };
	double rnorm{ L * std::pow(std::cos(lambda_final_rad), 2) / Re };
	double Br{ -B0 / std::pow(rnorm, 3) * 2 * sin(lambda_final_rad) };
	double Blambda{ -B0 / std::pow(rnorm, 3) * cos(lambda_final_rad) };

	return Br * std::cos(lambda_final_rad) + Blambda * std::sin(lambda_final_rad);
}

DLLEXPORT double getBy(double s, double dipoleEquatorialDist, double ILATdegrees)
{
	if (firstrun) { std::cout << "getBy: " << s << "  " << dipoleEquatorialDist << "  " << ILATdegrees << std::endl; firstrun = false; }
	double L{ Re / std::pow(std::cos(ILATdegrees * PI / 180), 2) };
	double lambda_final{ getLambdaAtS(s, L, ILATdegrees) };
	double lambda_final_rad{ lambda_final * PI / 180 };
	double rnorm{ L * std::pow(std::cos(lambda_final_rad), 2) / Re };
	double Br{ -B0 / std::pow(rnorm, 3) * 2 * sin(lambda_final_rad) };
	double Blambda{ -B0 / std::pow(rnorm, 3) * cos(lambda_final_rad) };

	return Br * std::sin(lambda_final_rad) - Blambda * std::cos(lambda_final_rad);
}

DLLEXPORT double Bx()
{

}

DLLEXPORT double By()
{

}

double foRungeKuttaBx(double* funcArg, int arrayLen, double** LUT, bool qsps, bool alfven)
{	// funcArg requirements: [t_RK = 0, y_0, ...] where t_RK = {0, h/2, h}, initial t_RK should be 0, this func will take care of the rest
	// dy / dt = f(t, y), y(t_0) = y_0
	// remaining funcArg elements are whatever you need in your callback function passed in
	//args array: [t_RKiter, vz, mu, q, m, pz_0, simtime, dt, omega E, const E]
	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];

	k1 = Bx(); //k1 = f(t_n, y_n), units of dy / dt

	funcArg[0] = funcArg[7] / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = Bx(); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = Bx(); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = funcArg[7];
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = Bx(); //k4 = f(t_n + h, y_n + h k3)

	return (k1 + 2 * k2 + 2 * k3 + k4) * funcArg[7] / 6; //returns delta y, not dy / dt, not total y
}

double foRungeKuttaBy(double* funcArg, int arrayLen, double** LUT, bool qsps, bool alfven)
{	// funcArg requirements: [t_RK = 0, y_0, ...] where t_RK = {0, h/2, h}, initial t_RK should be 0, this func will take care of the rest
	// dy / dt = f(t, y), y(t_0) = y_0
	// remaining funcArg elements are whatever you need in your callback function passed in
	//args array: [t_RKiter, vz, mu, q, m, pz_0, simtime, dt, omega E, const E]
	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];

	k1 = By(); //k1 = f(t_n, y_n), units of dy / dt

	funcArg[0] = funcArg[7] / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = By(); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = By(); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = funcArg[7];
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = By(); //k4 = f(t_n + h, y_n + h k3)

	return (k1 + 2 * k2 + 2 * k3 + k4) * funcArg[7] / 6; //returns delta y, not dy / dt, not total y
}

int main()
{
	double sIn{ 0.0 };
	double ILAT{ 72.0 };

	std::cout << "Specify s in meters: ";
	std::cin >> sIn;

	std::cout << "L: " << 1 / std::pow(std::cos(ILAT * PI / 180), 2) << std::endl;
	
	double L{ Re / std::pow(std::cos(ILAT * PI / 180), 2) };
	double lambda_final{ getLambdaAtS(sIn, L, ILAT) };

	std::cout << "\n================\nAngle Theta of location " << sIn << " meters along field line of ILAT " << ILAT << ": ";
	std::cout << std::fixed << std::setprecision(10) << lambda_final << "\n================" << std::endl;

	double lambda_final_rad{ lambda_final * PI / 180 };
	double r{ L * std::pow(std::cos(lambda_final_rad), 2) / Re };
	double Br{ -B0 / std::pow(r, 3) * 2 * sin(lambda_final_rad) };
	double Blambda{ -B0 / std::pow(r, 3) * cos(lambda_final_rad) };
	double B{ B0 / std::pow(r, 3) * std::sqrt(1 + 3 * std::pow(std::sin(lambda_final_rad), 2)) };
	double Balt{ std::sqrt(Br * Br + Blambda * Blambda) };

	std::cout << "lambda: " << lambda_final << ", r: " << r << ", Br: " << Br << ", Blambda: " << Blambda << ", B: " << B*1e9 << " nT, Balt: " << Balt << std::endl;

	return 0;
}