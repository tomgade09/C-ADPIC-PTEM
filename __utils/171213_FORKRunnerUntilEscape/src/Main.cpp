#include <cmath>
#include <iostream>

constexpr double MASS_PROTON{ 1.6726219e-27 }; //kg
constexpr double MASS_ELECTRON{ 9.1093836e-31 }; //kg
constexpr double CHARGE_ELEM{ 1.6021766e-19 }; //C
constexpr double RADIUS_EARTH{ 6.371e6 };		 //m
constexpr double BFIELD_EARTH{ -32.5e-6 };		 //T (at surface - 1 Re, Wiki mentioned a range from 25-65 uT, B0 would be about this, negative so B points into the Earth at North Pole)
constexpr double PI{ 3.1415927 };
constexpr double B0ATTHETA{ BFIELD_EARTH *  1.9102530 };
constexpr double DT{ 0.01 };
constexpr double MIN_Z_SIM{ (2.0e6 + RADIUS_EARTH) };
constexpr double MAX_Z_SIM{ 4 * RADIUS_EARTH };
constexpr double MIN_Z_NORM{ MIN_Z_SIM / RADIUS_EARTH };
constexpr double MAX_Z_NORM{ MAX_Z_SIM / RADIUS_EARTH };

double BFieldatZ(double z, double simtime)
{//for now, a simple dipole field
	if (z == 0)
		return 0.0; //add an error here if this case is true, at some point

	double norm{ RADIUS_EARTH };

	if ((z < RADIUS_EARTH) && (z > MIN_Z_NORM))
		norm = 1.0;

	return B0ATTHETA / pow(z / norm, 3); //Bz = B0 at theta * (1/rz(in Re))^3
}

void printArgs(double* args, int len)
{
	for (int iii = 0; iii < len; iii++)
		std::cout << args[iii] << ", ";
}

double accel1dCUDA(double* args, int len, double** LUT) //made to pass into 1D Fourth Order Runge Kutta code
{//args array: [t_RK, vz, mu, q, m, pz_0, simtime]
	//std::cout << "accel1dCUDA Array: ";
	//printArgs(args, len);
	//std::cout << "\n";

	double F_mir, ztmp;
	ztmp = args[5] + args[1] * args[0]; //pz_0 + vz * t_RK
	
	//Mirror force
	F_mir = -args[2] * B0ATTHETA * (-3 / (pow(ztmp / RADIUS_EARTH, 4)) / RADIUS_EARTH); //mu in [kg.m^2 / s^2.T] = [N.m / T]

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

	k1 = accel1dCUDA(funcArg, arrayLen, LUT); //k1 = f(t_n, y_n), units of dy / dt
	
	funcArg[0] = DT / 2;
	funcArg[1] = y_0 + k1 * funcArg[0];
	k2 = accel1dCUDA(funcArg, arrayLen, LUT); //k2 = f(t_n + h/2, y_n + h/2 k1)

	funcArg[1] = y_0 + k2 * funcArg[0];
	k3 = accel1dCUDA(funcArg, arrayLen, LUT); //k3 = f(t_n + h/2, y_n + h/2 k2)

	funcArg[0] = DT;
	funcArg[1] = y_0 + k3 * funcArg[0];
	k4 = accel1dCUDA(funcArg, arrayLen, LUT); //k4 = f(t_n + h, y_n + h k3)

	return (k1 + 2 * k2 + 2 * k3 + k4) * DT / 6; //returns delta y, not dy / dt, not total y
}

int main()
{
	double v_para{ 0.0 };
	double v_perp{ 0.0 };
	double z{ 0.0 };
	double simtime{ 0.0 };
	bool elecTF{ 1 };
	
	//accept input

	std::cout << "Input v_para (m/s): ";
	std::cin >> v_para;
	std::cout << "\nInput v_perp (m/s): ";
	std::cin >> v_perp;
	std::cout << "\nInput z (m): ";
	std::cin >> z;
	//std::cout << "\nInput electron?  >0 = true: ";
	//std::cin >> elecTF;

	std::cout << "Simulation Starting:\n" << "v_para: " << v_para << "\nv_perp: " << v_perp << "\nz: " << z << "\nsimtime: " << simtime << "\nelecTF: " << elecTF << "\n\n\n";

	//convert to mu
	v_perp = 0.5 * ((elecTF) ? (MASS_ELECTRON) : (MASS_PROTON)) * v_perp * v_perp / BFieldatZ(z, simtime);
	std::cout << "mu: " << v_perp << "\n\n\n";

	double args[7];
	args[3] = CHARGE_ELEM * ((elecTF) ? (-1.0) : (1.0));
	args[4] = (elecTF) ? (MASS_ELECTRON) : (MASS_PROTON);

	while ((z < MAX_Z_SIM * 1.001) && (z > MIN_Z_SIM * 0.999) && simtime <= 5000)
	{
		//[t_RK, vz, mu, q, m, pz_0, simtime]
		args[0] = 0.0;
		args[1] = v_para;
		args[2] = v_perp;
		args[5] = z;
		args[6] = simtime;
		
		double dv{ foRungeKuttaCUDA(args, 7, nullptr) };

		if (simtime == 0) { std::cout << dv << "    " << dv / DT << "    " << MASS_ELECTRON * dv / DT << "\n\n"; }

		v_para += dv;
		z += v_para * DT;
		simtime += DT;

		if (static_cast<int>(simtime/DT) % 10000000 == 0)
			std::cout << "Iterating: " << simtime << " s / 5000 s\n";
	}

	//convert to v_perp
	v_perp = sqrt(2 * v_perp * BFieldatZ(z, simtime) / args[4]);

	v_para /= 6.371e6;
	v_perp /= 6.371e6;
	z /= 6.371e6;

	if (simtime == 5000) { std::cout << "Didn't escape!! Time: " << simtime - DT << "\npara: " << v_para << "\nperp: " << v_perp << "\nz: " << z; }
	else { std::cout << "Escaped at time: " << simtime - DT << "\npara: " << v_para << "\nperp: " << v_perp << "\nz: " << z; }
	
	return 0;
}