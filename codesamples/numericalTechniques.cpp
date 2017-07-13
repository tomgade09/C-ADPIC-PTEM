#include <iostream>
#include "bfield\WireCoil_c.h"
#include <Eigen/Dense>
#include <cmath>

typedef Eigen::Vector3d (* FncPnt3dVector_t)(double*, int);

Eigen::Vector3d fourthOrderRungeKutta3DVector(FncPnt3dVector_t funcPointer, double* funcArg, int arrayLen, double h)
{	// funcArg requirements: [dt, var1, var2, var3...] where dt = {0, h/2, h}
	// var1, var2, var3 are the three components of Eigen::Vector3d inital values and values returned to tmp1, tmp2, tmp3, tmp4
	// dy / dt = f(t, y), y(t_0) = y_0 //(y_0 is initial [var1, var2, var3] -> funcArg[1], funcArg[2], funcArg[3])
	Eigen::Vector3d tmp1, tmp2, tmp3, tmp4, v_0;
	v_0 << funcArg[1], 
		   funcArg[2],
		   funcArg[3];

	tmp1 = funcPointer(funcArg, arrayLen); //k1 = f(t_n, y_n), units of dy / dt
	
	funcArg[0] = h / 2;
	funcArg[1] = v_0[0] + tmp1[0] * funcArg[0] / funcArg[5]; //not a big fan of this solution - cb function needs to return accel vs force
	funcArg[2] = v_0[1] + tmp1[1] * funcArg[0] / funcArg[5];
	funcArg[3] = v_0[2] + tmp1[2] * funcArg[0] / funcArg[5];
	tmp2 = funcPointer(funcArg, arrayLen); //k2 = f(t_n + h/2, y_n + h/2 k1)
	
	funcArg[1] = v_0[0] + tmp2[0] * funcArg[0] / funcArg[5];
	funcArg[2] = v_0[1] + tmp2[1] * funcArg[0] / funcArg[5];
	funcArg[3] = v_0[2] + tmp2[2] * funcArg[0] / funcArg[5];
	tmp3 = funcPointer(funcArg, arrayLen); //k3 = f(t_n + h/2, y_n + h/2 k2)
	
	funcArg[0] = h;
	funcArg[1] = v_0[0] + tmp3[0] * funcArg[0] / funcArg[5];
	funcArg[2] = v_0[1] + tmp3[1] * funcArg[0] / funcArg[5];
	funcArg[3] = v_0[2] + tmp3[2] * funcArg[0] / funcArg[5];
	tmp4 = funcPointer(funcArg, arrayLen); //k4 = f(t_n + h, y_n + h k3)
	
	return (tmp1 + 2 * tmp2 + 2 * tmp3 + tmp4) * h / 6; //returns units of y, not dy / dt
}

Eigen::Vector3d calcLorentz(Eigen::Vector3d p, Eigen::Vector3d v, double q, double m);
Eigen::Vector3d calcMirror(Eigen::Vector3d p, Eigen::Vector3d v, double q, double m, double scaledLength);
Eigen::Vector3d calcBatP(Eigen::Vector3d p);
Eigen::Vector3d calcEatP(Eigen::Vector3d p);

/*auto indexLambdaFcn = [=](int xsteps, int ysteps, int zsteps) -> unsigned int {
	return 3 * (zsteps * xsize * ysize + ysteps * xsize + xsteps); };*/

Eigen::Vector3d totalAccelCB(double* args, int len)
{//args array: [dt, vx, vy, vz, q, m, px_0, py_0, pz_0]
	if (len != 9)
	{
		std::cout << "Array is not the right length.  Proper array format is: [dt, vx, vy, vz, q, m, px_0, py_0, pz_0].  Returning zero vector.\n";
		Eigen::Vector3d ret;
		ret << 0., 0., 0.;
		return ret;
	}
	
	Eigen::Vector3d F_lor, F_mir, p, v;
	p << args[6] + args[1] * args[0], 
		 args[7] + args[2] * args[0],
		 args[8] + args[3] * args[0];
	v << args[1], 
		 args[2],
		 args[3];
	F_lor = calcLorentz(p, v, args[4], args[5]);
	F_mir = calcMirror(p, v, args[4], args[5], 1e-8); //Not sure if this is a good scaled length for these sims - prob not
	
	return (F_lor + F_mir) / args[5];
}//calculate force: first calc B (for Lorentz), then slope B (for mirror), needs to return accel

Eigen::Vector3d calcLorentz(Eigen::Vector3d p, Eigen::Vector3d v, double q, double m)
{
	Eigen::Vector3d B, E;
	B = calcBatP(p);
	E = calcEatP(p);
	
	return (q * (E + v.cross(B)));
}

Eigen::Vector3d calcMirror(Eigen::Vector3d p, Eigen::Vector3d v, double q, double m, double scaledLength)
{
	Eigen::Vector3d Bp, B_plus, B_minus, halfds, F_gradB, p_minus, p_plus;
	
	Bp = calcBatP(p);
	double Bp_len{ std::sqrt(Bp.dot(Bp)) };
	halfds = Bp * scaledLength / Bp_len;
	
	p_minus = p - halfds;
	p_plus = p + halfds;
	B_minus = calcBatP(p - halfds);
	B_plus = calcBatP(p + halfds);

	double vperp2 = v.dot(v) - pow(Bp.dot(v) / Bp_len, 2);
	double mu = m * vperp2 / (2 * Bp_len);

	F_gradB = (B_plus - B_minus) * - mu;

	for (unsigned int jjj = 0; jjj < 3; jjj++)
	{
		if (abs(Bp[jjj]) < 1e-20) //do I still need this "correction"?  As is, it won't even run (B will be bigger than 1e-20)
		{
			F_gradB[jjj] = 0;
			continue;
		}
	}
	
	F_gradB /= (2 * scaledLength);

	return F_gradB;
}