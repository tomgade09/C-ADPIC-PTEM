#include <cmath>
#include <iostream>
#include <vector>
#include "FileIO\fileIO.h"

constexpr double MASS_PROTON{ 1.6726219e-27 }; //kg
constexpr double MASS_ELECTRON{ 9.1093836e-31 }; //kg
constexpr double CHARGE_ELEM{ 1.6021766e-19 }; //C
constexpr double RADIUS_EARTH{ 6.371e6 };		 //m
constexpr double BFIELD_EARTH{ -32.5e-6 };		 //T (at surface - 1 Re, Wiki mentioned a range from 25-65 uT, B0 would be about this, negative so B points into the Earth at North Pole)
constexpr double PI{ 3.1415927 };
constexpr double B0ATTHETA{ BFIELD_EARTH *  1.9102530 };
constexpr double GLOBAL_DT{ 0.1 };
constexpr double MIN_Z_SIM{ (2.0e6 + RADIUS_EARTH) };
constexpr double MAX_Z_SIM{ 4 * RADIUS_EARTH };
constexpr double MIN_Z_NORM{ MIN_Z_SIM / RADIUS_EARTH };
constexpr double MAX_Z_NORM{ MAX_Z_SIM / RADIUS_EARTH };
constexpr double ERR_LOW_THRESHOLD{ 0.0000001 }; //allow for up to a 1 % error
constexpr bool elecTF{ 1 };
double DT{ 0.1 };

double BFieldatZ(double z, double simtime)
{//for now, a simple dipole field
	if (z == 0)
		return 0.0; //add an error here if this case is true, at some point

	double norm{ RADIUS_EARTH };

	if ((z < RADIUS_EARTH / 1e5) && (z > MIN_Z_NORM))
		norm = 1.0;

	return B0ATTHETA / pow(z / norm, 3); //Bz = B0 at theta * (1/rz(in Re))^3
}

double gradB(double z, double simtime)
{
	if (z == 0)
		return 0.0;

	double norm{ RADIUS_EARTH };

	if ((z < RADIUS_EARTH / 1e5) && (z > MIN_Z_NORM))
		norm = 1.0;

	return B0ATTHETA * (-3 / (pow(z / norm, 4))) / RADIUS_EARTH;
}

void printArgs(double* args, int len)
{
	for (int iii = 0; iii < len; iii++)
		std::cout << args[iii] << ", ";
}

double accel1dCUDA(double* args, int len, double** LUT) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [t_RK, vz, mu, q, m, pz_0, simtime]
	//std::cout << "\naccel1dCUDA Array: ";
	//printArgs(args, len);
	//std::cout << "\n";

	double F_mir, ztmp;
	ztmp = args[5] + args[1] * args[0]; //pz_0 + vz * t_RK
	//std::cout << "ztmp: " << ztmp;
	
	//Mirror force
	F_mir = -args[2] * gradB(ztmp, args[6]); //mu in [kg.m^2 / s^2.T] = [N.m / T]
	//std::cout << " Grad B: " << gradB(ztmp, args[6]) << "  F_mir: " << F_mir << " / " << args[4] << "\n\n";

	return (F_mir) / args[4];
}//returns an acceleration in the parallel direction to the B Field

double foRungeKuttaCUDA(double* funcArg, int arrayLen, double** LUT)
{	// funcArg requirements: [t_RK = 0, y_0, ...] where t_RK = {0, h/2, h}, initial t_RK should be 0, this func will take care of the rest
	// dy / dt = f(t, y), y(t_0) = y_0
	// remaining funcArg elements are whatever you need in your callback function passed in
	//std::cout << "foRungeKuttaCUDA Array: ";
	//printArgs(funcArg, arrayLen);
	//std::cout << "\n";

	double k1, k2, k3, k4, y_0;
	y_0 = funcArg[1];
	//std::cout << "y_0: " << y_0 << "\n";
	k1 = accel1dCUDA(funcArg, arrayLen, LUT); //k1 = f(t_n, y_n), units of dy / dt
	//std::cout << "k1: " << k1 << "\n";
	
	//std::cout << "t_RK from " << funcArg[0];
	funcArg[0] = DT / 2;
	//std::cout << " to " << funcArg[0] << " vz from " << funcArg[1] << " to ";
	funcArg[1] = y_0 + k1 * funcArg[0];
	//std::cout << funcArg[1] << "\n";
	k2 = accel1dCUDA(funcArg, arrayLen, LUT); //k2 = f(t_n + h/2, y_n + h/2 k1)
	//std::cout << "k2: " << k2 << "\n";
	
	//std::cout << "t_RK " << funcArg[0] << " vz from " << funcArg[1] << " to ";
	funcArg[1] = y_0 + k2 * funcArg[0];
	//std::cout << funcArg[1] << "\n";
	k3 = accel1dCUDA(funcArg, arrayLen, LUT); //k3 = f(t_n + h/2, y_n + h/2 k2)
	//std::cout << "k3: " << k3 << "\n";
	
	//std::cout << "t_RK from " << funcArg[0];
	funcArg[0] = DT;
	//std::cout << " to " << funcArg[0] << " vz from " << funcArg[1] << " to ";
	funcArg[1] = y_0 + k3 * funcArg[0];
	//std::cout << funcArg[1] << "\n";
	k4 = accel1dCUDA(funcArg, arrayLen, LUT); //k4 = f(t_n + h, y_n + h k3)
	//std::cout << "k4: " << k4 << "\n";

	//std::cout << "Complete with a FORK run\nAcceleration: " << (k1 + 2 * k2 + 2 * k3 + k4) / 6 << "\nv_tot: " << y_0 + (k1 + 2 * k2 + 2 * k3 + k4) * DT / 6 << "\n\n";

	return (k1 + 2 * k2 + 2 * k3 + k4) * DT / 6; //returns delta y, not dy / dt, not total y
}

double** runFORK(double v_para, double mu, double z, double timetrigger, int power)
{
	//std::cout << "Write initial values to vector   " << DT << "   " << static_cast<int>(timetrigger * pow(10, power)) + 1 << "\n";
	//std::vector<std::vector<double>> attributes(4, std::vector<double>((static_cast<int>(timetrigger * pow(10, power)) + 1), 0));

	double** attributes = new double*[4];
	//std::cout << 21 << "  " << timetrigger << "\n";
	for (int iii = 0; iii < 4; iii++)
		attributes[iii] = new double[21];

	double simtime{ 0.0 };
	long long loopind{ 0 };
	int ind{ 0 };

	long long dtsPerGlobalDT{ static_cast<long long>(GLOBAL_DT / DT) };
	std::cout << dtsPerGlobalDT << "\n";

	double args[7];
	args[3] = CHARGE_ELEM * ((elecTF) ? (-1.0) : (1.0));
	args[4] = (elecTF) ? (MASS_ELECTRON) : (MASS_PROTON);
	
	attributes[0][ind] = simtime;
	attributes[1][ind] = v_para;
	attributes[2][ind] = sqrt(2 * mu * BFieldatZ(z, simtime) / args[4]);
	attributes[3][ind] = z;
	
	while ((simtime < timetrigger)) //&& (z < (MAX_Z_SIM * 1.0001)) && (z > (MIN_Z_SIM * 0.9999)))
	{
		//std::cout << "iterate " << simtime <<"\n";
		//[t_RK, vz, mu, q, m, pz_0, simtime]
		args[0] = 0.0;
		args[1] = v_para;
		args[2] = mu;
		args[5] = z;
		args[6] = 0.0;

		//std::cout << "[t_RK,  vz,  mu,  q,  m,  pz_0,  simtime]\n";

		v_para += foRungeKuttaCUDA(args, 7, nullptr);
		z += v_para * DT;
		simtime += DT;
		loopind++;

		if (loopind % dtsPerGlobalDT == 0)
		{
			//std::cout << simtime << "\n";
			ind++;

			if (ind > 21)
			{
				std::cout << ind << "  " << simtime << "  Overflowed buffer.  Returning.  Results may be suspect.\n";
				return attributes;
			}

			attributes[0][ind] = simtime;
			attributes[1][ind] = v_para;
			attributes[2][ind] = sqrt(2 * mu * BFieldatZ(z, simtime) / args[4]);
			attributes[3][ind] = z;
		}
	}

	std::cout << "DT:    " << DT << "    Complete at " << simtime <<"!\n";
	return attributes;
}

int main()
{
	double v_para{ 0.0 };
	double v_perp{ 0.0 };
	double z{ 0.0 };

	//std::vector<std::vector<std::vector<double>>> results;
	//results.reserve(10);
	double*** results = new double**[10];

	//accept input
	std::cout << "Input v_para (m/s): ";
	std::cin >> v_para;
	std::cout << "\nInput v_perp (m/s): ";
	std::cin >> v_perp;
	std::cout << "\nInput z (m): ";
	std::cin >> z;

	std::cout << "Simulation Starting:\n" << "v_para: " << v_para << "\nv_perp: " << v_perp << "\nz: " << z << "\nelecTF: " << elecTF << "\n";

	//convert to mu
	v_perp = 0.5 * ((elecTF) ? (MASS_ELECTRON) : (MASS_PROTON)) * v_perp * v_perp / BFieldatZ(z, 0.0);
	std::cout << "mu: " << v_perp << "\n\n\n";

	int iii{ 0 };
	int power{ -1 };
	bool condition{ true };
	while (condition)
	{
		DT = pow(10, power--);
		std::cout << "Start DT: " << DT << "\n";
		//results.push_back(runFORK(v_para, v_perp, z, 0.02, -power));
		results[iii] = runFORK(v_para, v_perp, z, 2.0, power);

		//long highDTind  { static_cast<long>(results.size() - 1) }; //index of the latest DT FORK data - ex: second iteration this value is be 1 - this data is at results[1]
		//long highMsmtInd{ static_cast<long>(results[highDTind][0].size() - 1) }; //index of the latest data's last element - ex: first measurement (0.01s) with a timecap of 0.02s - this value is 2 - data is at results[0][attr][2]
		//long nextMsmtInd{ (highDTind > 0) ? static_cast<long>(results[highDTind - 1][0].size() - 1) : 0 }; //index of one previous data's last element - ex: two measurements have taken place, this refers to the first and is (results[0][0].size - 1)
		//std::cout << highDTind << "  " << highMsmtInd << "  " << nextMsmtInd << "\n";

		//std::cout << "condition check.  Hold onto your butts...\n";
		//std::cout << "true!\n";
		//int highind{ static_cast<int>(2 * pow(10, iii + 2)) };
		//int lastind{ static_cast<int>(2 * pow(10, iii + 1)) };
		int highind{ 20 };
		int lastind{ 20 };


		if (iii > 0)
		{
			double vrecent{ results[iii][1][highind] };
			double vlast{ results[iii - 1][1][lastind] };
			std::cout << vlast << " / " << vrecent << " = " << abs(1 - vlast / vrecent) << "\n";
			//std::cout << "power: " << power << "\n";
			if (power < -6)//((abs(1 - vlast / vrecent) < ERR_LOW_THRESHOLD) || (power < -8))
			{
				std::cout << "condition = true\n";
				condition = false;
			}
		}

		iii++;
	}

	std::cout << "Done with FORK loop!\n";

	for (int jjj = 0; jjj < iii; jjj++) //measurements
	{
		double* data[4]{ results[jjj][0], results[jjj][1], results[jjj][2], results[jjj][3] };
		fileIO::write2DCSV(("./../../../out/1e" + std::to_string(-1 - jjj) + ".csv").c_str(), data, 21, 4, ',');
	}

	return 0;
}