#include <iostream>
#include "bfield/WireCoil_c.h"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include "include/numericaltechniques.h"

//needs to be templated and go into a header
Eigen::Vector3d fourthOrderRungeKutta3DVector(FncPnt3dVector_t funcPointer, double* funcArg, int arrayLen, double h)
{	// funcArg requirements: [dt, var1, var2, var3...] where dt = {0, h/2, h}, initial dt should be 0, this func will take care of the rest
	// var1, var2, var3 are the three components of Eigen::Vector3d inital values and values returned to k1, k2, k3, k4
	// dy / dt = f(t, y), y(t_0) = y_0 //(y_0 is initial [var1, var2, var3] -> funcArg[1], funcArg[2], funcArg[3])
	// remaining funcArg elements are whatever you need in your callback function passed in
	Eigen::Vector3d k1, k2, k3, k4, v_0;
	v_0 << funcArg[1], 
		   funcArg[2],
		   funcArg[3];

	k1 = funcPointer(funcArg, arrayLen); //k1 = f(t_n, y_n), units of dy / dt
	
	funcArg[0] = h / 2;
	funcArg[1] = v_0[0] + k1[0] * funcArg[0];
	funcArg[2] = v_0[1] + k1[1] * funcArg[0];
	funcArg[3] = v_0[2] + k1[2] * funcArg[0];
	k2 = funcPointer(funcArg, arrayLen); //k2 = f(t_n + h/2, y_n + h/2 k1)
	
	funcArg[1] = v_0[0] + k2[0] * funcArg[0];
	funcArg[2] = v_0[1] + k2[1] * funcArg[0];
	funcArg[3] = v_0[2] + k2[2] * funcArg[0];
	k3 = funcPointer(funcArg, arrayLen); //k3 = f(t_n + h/2, y_n + h/2 k2)
	
	funcArg[0] = h;
	funcArg[1] = v_0[0] + k3[0] * funcArg[0];
	funcArg[2] = v_0[1] + k3[1] * funcArg[0];
	funcArg[3] = v_0[2] + k3[2] * funcArg[0];
	k4 = funcPointer(funcArg, arrayLen); //k4 = f(t_n + h, y_n + h k3)
	
	return (k1 + 2 * k2 + 2 * k3 + k4) * h / 6; //returns units of y, not dy / dt
}

//this one can't really be templated because of p and v
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
	F_lor = calcLorentz(p, v, args[4]);
	F_mir = calcMirror(p, v, args[4], args[5], 1e-8); //Not sure if this is a good scaled length for these sims - prob not -> total z range is many km
	//F_mir << 0., 0., 0.;

	return (F_lor + F_mir) / args[5];
}//calculate force: first calc B (for Lorentz), then slope B (for mirror), needs to return accel

//template this function - need a 2D case, not sure how to cross product
Eigen::Vector3d calcLorentz(Eigen::Vector3d p, Eigen::Vector3d v, double q)
{//returns force - needs a wrapper function to divide by m: 4 order RK [k_n] = [dy / dt], in this case [k_n] = [acceleration]
	Eigen::Vector3d B, E;
	B = calcBatP(p);
	E = calcEatP(p);
	
	return (q * (E + v.cross(B)));
}

//template this function
Eigen::Vector3d calcMirror(Eigen::Vector3d p, Eigen::Vector3d v, double q, double m, double scaledLength)
{//returns force - needs a wrapper function to divide by m before going back to 4RK: 4 order RK [k_n] = [dy / dt], in this case [k_n] = [acceleration]
	Eigen::Vector3d Bp, B_plus, B_minus, halfds, F_gradB;// , p_minus, p_plus;
	
	Bp = calcBatP(p);
	double Bp_len{ std::sqrt(Bp.dot(Bp)) };
	halfds = Bp * scaledLength / Bp_len;
	
	B_minus = calcBatP(p - halfds);
	B_plus = calcBatP(p + halfds);

	double vperp2 = v.dot(v) - pow(Bp.dot(v) / Bp_len, 2);
	double mu = m * vperp2 / (2 * Bp_len);

	F_gradB = (B_plus - B_minus) * - mu;

	/*for (unsigned int jjj = 0; jjj < 3; jjj++) //not general enough
	{
		if (abs(Bp[jjj]) < 1e-20) //do I still need this "correction"?  As is, it won't even run (B will be bigger than 1e-20).
		{						  //needs to be removed/updated due to above, but I'll leave for now
			F_gradB[jjj] = 0;
			continue;
		}
	}*/
	
	F_gradB /= (2 * scaledLength);

	return F_gradB;
}