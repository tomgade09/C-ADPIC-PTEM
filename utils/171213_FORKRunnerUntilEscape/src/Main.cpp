#include <cmath>
#include <iostream>

constexpr double MASS_PROTON{ 1.672621898e-27 };	//kg
constexpr double MASS_ELECTRON{ 9.10938356e-31 };	//kg
constexpr double CHARGE_ELEM{ 1.6021766209e-19 }; //C
constexpr double RADIUS_EARTH{ 6.3712e6 };			//m
constexpr double JOULE_PER_EV{ 1.6021766209e-19 };
constexpr double PI{ 3.14159265358979323846 };

constexpr double DT{ 0.01 };
constexpr double MIN_Z_SIM{ 2030837.49610366 };
constexpr double MAX_Z_SIM{ 19881647.2473464 };
constexpr double MIN_Z_NORM{ MIN_Z_SIM / RADIUS_EARTH };
constexpr double MAX_Z_NORM{ MAX_Z_SIM / RADIUS_EARTH };

double sinpi(double x)
{
	return sin(PI * x);
}

double cospi(double x)
{
	return cos(PI * x);
}

//B Field related kernels
double getSAtLambda(double* consts, int arrayLength, double lambdaDegrees)///FIX TO GIVE S FROM RE NOT EQUATOR!!!!!!!!!!!AA!!!1111!1!!111!
{
	double x{ asinh(sqrt(3.0) * sinpi(lambdaDegrees / 180.0)) };

	return (0.5 * consts[2] / sqrt(3.0)) * (x + sinh(x) * cosh(x)); /* L */
}

double getLambdaAtS(double* consts, int arrayLength, double s)
{// consts: [ B0, ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	double lambda_tmp{ (-consts[1] / consts[4]) * s + consts[1] }; //-ILAT / s_max * s + ILAT
	double s_tmp{ consts[4] - getSAtLambda(consts, arrayLength, lambda_tmp) };
	double dlambda{ 1.0 };
	bool   over{ 0 };

	while (abs((s_tmp - s) / s) > consts[6]) //errorTolerance
	{
		while (1)
		{
			over = (s_tmp >= s);
			if (over)
			{
				lambda_tmp += dlambda;
				s_tmp = consts[4] - getSAtLambda(consts, arrayLength, lambda_tmp);
				if (s_tmp < s)
					break;
			}
			else
			{
				lambda_tmp -= dlambda;
				s_tmp = consts[4] - getSAtLambda(consts, arrayLength, lambda_tmp);
				if (s_tmp >= s)
					break;
			}
		}
		if (dlambda < consts[6] / 100.0) //errorTolerance
			break;
		dlambda /= 5.0; //through trial and error, this reduces the number of calculations usually (compared with 2, 2.5, 3, 4, 10)
	}

	return lambda_tmp;
}

double BFieldAtS(double* consts, int arrayLength, double s, double simtime)
{// consts: [ B0, ILATDeg, L, L_norm, s_max, ds, errorTolerance ]
	double lambda_deg{ getLambdaAtS(consts, arrayLength, s) };
	double rnorm{ consts[3] * cospi(lambda_deg / 180.0) * cospi(lambda_deg / 180.0) };

	return -consts[0] / (rnorm * rnorm * rnorm) * sqrt(1.0 + 3 * sinpi(lambda_deg / 180.0) * sinpi(lambda_deg / 180.0));
}

double gradBAtS(double* consts, int arrayLength, double s, double simtime)
{
	return (BFieldAtS(consts, arrayLength, s + consts[5], simtime) - BFieldAtS(consts, arrayLength, s - consts[5], simtime)) / (2 * consts[5]);
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

	double Bargs[]{ 3.12e-5, 72.0, RADIUS_EARTH / pow(cospi(72.0 / 180.0), 2), 1 / pow(cospi(72.0 / 180.0), 2), 0.0, 6371.2, 10e-4 };
	Bargs[4] = getSAtLambda(Bargs, 7, 72.0);

	double F_mir, ztmp;
	ztmp = args[5] + args[1] * args[0]; //pz_0 + vz * t_RK
	
	//Mirror force
	F_mir = -args[2] * gradBAtS(Bargs, 7, ztmp, 0.0); //mu in [kg.m^2 / s^2.T] = [N.m / T]

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

	double Bargs[]{ 3.12e-5, 72.0, RADIUS_EARTH / pow(cospi(72.0 / 180.0), 2), 1 / pow(cospi(72.0 / 180.0), 2), 0.0, 6371.2, 10e-4 };
	Bargs[4] = getSAtLambda(Bargs, 7, 72.0);

	std::cout << "\n\nB Field At Top: " << BFieldAtS(Bargs, 7, MAX_Z_SIM, 0.0) << "\nB Field At Btm: " << BFieldAtS(Bargs, 7, MIN_Z_SIM, 0.0) << std::endl;

	//convert to mu
	v_perp = 0.5 * ((elecTF) ? (MASS_ELECTRON) : (MASS_PROTON)) * v_perp * v_perp / BFieldAtS(Bargs, 7, z, simtime);
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
	v_perp = sqrt(2 * v_perp * BFieldAtS(Bargs, 7, z, simtime) / args[4]);

	v_para /= RADIUS_EARTH;
	v_perp /= RADIUS_EARTH;
	z /= RADIUS_EARTH;

	if (simtime == 5000) { std::cout << "Didn't escape!! Time: " << simtime - DT << "\npara: " << v_para << "\nperp: " << v_perp << "\nz: " << z; }
	else { std::cout << "Escaped at time: " << simtime - DT << "\npara: " << v_para << "\nperp: " << v_perp << "\nz: " << z; }
	
	return 0;
}